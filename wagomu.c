/*
 * Copyright (C) 2009 The Tegaki project contributors
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/*
 * Contributors to this file:
 *  - Mathieu Blondel
 *  - Emanuel Ferraz
 */

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "wagomu.h"

#define MAGIC_NUMBER 0x77778888
#define VEC_DIM_MAX 4
#define NORM_SIZE 1000.0f
#define ARENA_DEFAULT_SIZE (1024 * 1024)

#undef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#undef MIN3
#define MIN3(a, b, c) (MIN((a), MIN((b), (c))))

/**
 * Type definitions.
 */

typedef struct {
  char *base;
  unsigned int size;
  unsigned int offset;
} wagomu_arena_t;

typedef struct {
  unsigned int unicode;
  unsigned int n_vectors;
} wagomu_character_info_t;

typedef struct {
  unsigned int n_strokes;
  unsigned int n_chars;
  unsigned int offset;
  char pad[4];
} wagomu_character_group_t;

typedef struct {
  float x;
  float y;
} wagomu_point_t;

typedef struct {
  unsigned int start;
  unsigned int count;
} wagomu_stroke_meta_t;

typedef struct {
  wagomu_point_t *points;
  unsigned int points_size;
  unsigned int points_cap;

  wagomu_stroke_meta_t *meta;
  unsigned int n_strokes;
  unsigned int meta_cap;
} wagomu_strokes_t;

struct wagomu_recognizer_s {
  wagomu_arena_t arena;

  const char *data;
  const wagomu_character_info_t *characters;
  const wagomu_character_group_t *groups;
  float *dtw1;
  float *dtw2;
  char *error_msg;
  wagomu_result_t *distm;
  wagomu_strokes_t strokes;
  wagomu_strokes_t proc_strokes;
  float *features;

  unsigned int features_cap;
  unsigned int n_characters;
  unsigned int n_groups;
  unsigned int dimension;
  unsigned int downsample_threshold;
  unsigned int window_size;
};

/*
 * Arena allocator functions.
 */

static wagomu_arena_t arena_create(unsigned int size) {
  wagomu_arena_t arena;
  arena.base = (char *)malloc(size);
  arena.size = size;
  arena.offset = 0;
  return arena;
}

static void *arena_alloc(wagomu_arena_t *arena, unsigned int size,
                         unsigned int align) {
  unsigned int aligned_offset = (arena->offset + align - 1) & ~(align - 1);

  if (aligned_offset + size > arena->size) {
    unsigned int new_size = arena->size * 2;
    while (aligned_offset + size > new_size)
      new_size *= 2;

    char *new_base = (char *)realloc(arena->base, new_size);
    if (!new_base)
      return NULL;

    arena->base = new_base;
    arena->size = new_size;
  }

  void *ptr = arena->base + aligned_offset;
  arena->offset = aligned_offset + size;
  return ptr;
}

static void arena_destroy(wagomu_arena_t *arena) {
  if (arena->base)
    free(arena->base);

  arena->base = NULL;
  arena->size = 0;
  arena->offset = 0;
}

/*
 * Helper functions.
 */

// `strokes_start_stroke` starts a new stroke in the strokes container.
static void strokes_start_stroke(wagomu_arena_t *arena,
                                 wagomu_strokes_t *strokes) {
  if (strokes->n_strokes == strokes->meta_cap) {
    unsigned int new_cap = (strokes->meta_cap == 0) ? 4 : strokes->meta_cap * 2;
    wagomu_stroke_meta_t *new_meta = (wagomu_stroke_meta_t *)arena_alloc(
        arena, new_cap * sizeof(wagomu_stroke_meta_t),
        _Alignof(wagomu_stroke_meta_t));

    if (strokes->meta && strokes->n_strokes > 0)
      memcpy(new_meta, strokes->meta,
             strokes->n_strokes * sizeof(wagomu_stroke_meta_t));

    strokes->meta = new_meta;
    strokes->meta_cap = new_cap;
  }

  strokes->meta[strokes->n_strokes].start = strokes->points_size;
  strokes->meta[strokes->n_strokes].count = 0;
  strokes->n_strokes++;
}

// `strokes_add_point` adds a point to the current stroke.
static void strokes_add_point(wagomu_arena_t *arena, wagomu_strokes_t *strokes,
                              wagomu_point_t point) {
  if (strokes->n_strokes == 0)
    return;

  if (strokes->points_size == strokes->points_cap) {
    unsigned int new_cap =
        (strokes->points_cap == 0) ? 64 : strokes->points_cap * 2;
    wagomu_point_t *new_points = (wagomu_point_t *)arena_alloc(
        arena, new_cap * sizeof(wagomu_point_t), _Alignof(wagomu_point_t));

    if (strokes->points && strokes->points_size > 0)
      memcpy(new_points, strokes->points,
             strokes->points_size * sizeof(wagomu_point_t));

    strokes->points = new_points;
    strokes->points_cap = new_cap;
  }

  strokes->points[strokes->points_size++] = point;
  strokes->meta[strokes->n_strokes - 1].count++;
}

// `strokes_copy` deep copies source strokes into destination.
static void strokes_copy(wagomu_arena_t *arena, wagomu_strokes_t *dest,
                         const wagomu_strokes_t *src) {
  if (src->points_size == 0) {
    dest->points_size = 0;
    dest->n_strokes = 0;
    return;
  }

  if (dest->points_cap < src->points_size) {
    dest->points = (wagomu_point_t *)arena_alloc(
        arena, src->points_size * sizeof(wagomu_point_t),
        _Alignof(wagomu_point_t));
    dest->points_cap = src->points_size;
  }

  if (dest->meta_cap < src->n_strokes) {
    dest->meta = (wagomu_stroke_meta_t *)arena_alloc(
        arena, src->n_strokes * sizeof(wagomu_stroke_meta_t),
        _Alignof(wagomu_stroke_meta_t));
    dest->meta_cap = src->n_strokes;
  }

  memcpy(dest->points, src->points, src->points_size * sizeof(wagomu_point_t));
  memcpy(dest->meta, src->meta, src->n_strokes * sizeof(wagomu_stroke_meta_t));
  dest->points_size = src->points_size;
  dest->n_strokes = src->n_strokes;
}

// `dist_sq` returns the euclidean distance squared between two points.
static float dist_sq(wagomu_point_t point1, wagomu_point_t point2) {
  return (point1.x - point2.x) * (point1.x - point2.x) +
         (point1.y - point2.y) * (point1.y - point2.y);
}

// `normalize_strokes` resizes the strokes into a box of size
// `NORM_SIZE`x`NORM_SIZE` and centers them.
static void normalize_strokes(wagomu_strokes_t *strokes) {
  if (strokes->points_size == 0)
    return;

  float min_x = FLT_MAX, max_x = -FLT_MAX;
  float min_y = FLT_MAX, max_y = -FLT_MAX;

  for (unsigned int i = 0; i < strokes->points_size; ++i) {
    if (strokes->points[i].x < min_x)
      min_x = strokes->points[i].x;
    if (strokes->points[i].x > max_x)
      max_x = strokes->points[i].x;
    if (strokes->points[i].y < min_y)
      min_y = strokes->points[i].y;
    if (strokes->points[i].y > max_y)
      max_y = strokes->points[i].y;
  }

  float width = max_x - min_x;
  float height = max_y - min_y;
  if (width == 0)
    width = 1;
  if (height == 0)
    height = 1;

  float ratio_w = NORM_SIZE / width;
  float ratio_h = NORM_SIZE / height;
  float ratio = (ratio_w < ratio_h) ? ratio_w : ratio_h;

  float new_width = width * ratio;
  float new_height = height * ratio;
  float dx = (NORM_SIZE - new_width) / 2.0f;
  float dy = (NORM_SIZE - new_height) / 2.0f;

  for (unsigned int i = 0; i < strokes->points_size; ++i) {
    strokes->points[i].x = (strokes->points[i].x - min_x) * ratio + dx;
    strokes->points[i].y = (strokes->points[i].y - min_y) * ratio + dy;
  }
}

// `downsample_strokes` reduces point count based on distance.
static void downsample_strokes(wagomu_strokes_t *strokes,
                               unsigned int downsample_threshold) {
  if (strokes->points_size == 0)
    return;

  unsigned int write_idx = 0;
  float threshold_sq = (float)(downsample_threshold * downsample_threshold);

  for (unsigned int s = 0; s < strokes->n_strokes; ++s) {
    unsigned int start = strokes->meta[s].start;
    unsigned int count = strokes->meta[s].count;

    if (count == 0)
      continue;

    unsigned int new_start = write_idx;
    unsigned int new_count = 0;

    strokes->points[write_idx++] = strokes->points[start];
    new_count++;

    wagomu_point_t last_original = strokes->points[start + count - 1];

    for (unsigned int j = 1; j < count; ++j) {
      wagomu_point_t last_kept = strokes->points[write_idx - 1];
      if (dist_sq(strokes->points[start + j], last_kept) >= threshold_sq) {
        strokes->points[write_idx++] = strokes->points[start + j];
        new_count++;
      }
    }

    wagomu_point_t last_kept = strokes->points[write_idx - 1];
    if (count > 1 && dist_sq(last_original, last_kept) > 100.0f) {
      strokes->points[write_idx++] = last_original;
      new_count++;
    }

    strokes->meta[s].start = new_start;
    strokes->meta[s].count = new_count;
  }

  strokes->points_size = write_idx;
}

// tbh I don't really know the algorithm to know what this does.
static inline float local_distance(wagomu_recognizer_t *recognizer, float *vec1,
                                   float *vec2) {
  float sum = 0;
  for (unsigned int i = 0; i < recognizer->dimension; i++)
    sum += fabsf(vec2[i] - vec1[i]);

  return sum;
}

// `dtw` is the Dynamic Time Warping algorithm to measure similarity between two
// temporal sequences.
static inline float dtw(wagomu_recognizer_t *recognizer, float *seq1,
                        unsigned int n, float *seq2, unsigned int m) {
  float *seq2_start = seq2;

  recognizer->dtw1[0] = 0;
  recognizer->dtw2[0] = FLT_MAX;

  /* Initialize the edge cells */
  for (unsigned int i = 1; i < m; i++)
    recognizer->dtw1[i] = FLT_MAX;

  seq1 += recognizer->dimension;

  /* Iterate over columns */
  for (unsigned int i = 1; i < n; i++) {
    seq2 = seq2_start + recognizer->dimension;

    /* Iterate over cells of that column */
    for (unsigned int j = 1; j < m; j++) {
      float cost = local_distance(recognizer, seq1, seq2);
      /* Inductive step */
      recognizer->dtw2[j] =
          cost + MIN3(recognizer->dtw2[j - 1], recognizer->dtw1[j],
                      recognizer->dtw1[j - 1]);

      seq2 += recognizer->dimension;
    }

    float *tmp = recognizer->dtw1;
    recognizer->dtw1 = recognizer->dtw2;
    recognizer->dtw2 = tmp;
    *recognizer->dtw2 = FLT_MAX;

    seq1 += recognizer->dimension;
  }

  return recognizer->dtw1[m - 1];
}

// `char_dist_cmp` is a callback function for glibc's `qsort` function.
static int char_dist_cmp(const void *a, const void *b) {
  const wagomu_result_t *da = (const wagomu_result_t *)a;
  const wagomu_result_t *db = (const wagomu_result_t *)b;
  if (da->dist < db->dist)
    return -1;
  if (da->dist > db->dist)
    return 1;
  return 0;
}

// recognizer implementation.

// `wagomu_recognizer_new` allocates a new `wagomu_recognizer_t` on the heap
// and initializes it with the model from the given `path`.
// It will always return a struct, but the struct may not be completely
// initialized.
// Always check `wagomu_get_error_message` to see if there was
// an error or not.
// `model_bytes` is expected to have a lifetime greater or equal than
// the recognizer.
wagomu_recognizer_t *wagomu_recognizer_new(const char *model_bytes,
                                           unsigned int model_size) {
  wagomu_arena_t arena = arena_create(ARENA_DEFAULT_SIZE);
  if (!arena.base)
    return NULL;

  wagomu_recognizer_t *recognizer = (wagomu_recognizer_t *)arena_alloc(
      &arena, sizeof(wagomu_recognizer_t), _Alignof(wagomu_recognizer_t));

  memset(recognizer, 0, sizeof(wagomu_recognizer_t));
  recognizer->arena = arena;

  if (!model_bytes) {
    recognizer->error_msg = "model bytes cannot be null";
    return recognizer;
  }

  if (model_size < sizeof(unsigned int) * 4) {
    recognizer->error_msg = "Model file does not contain a header";
    return recognizer;
  }

  recognizer->data = model_bytes;

  unsigned int *header = (unsigned int *)recognizer->data;

  if (header[0] != MAGIC_NUMBER) {
    recognizer->error_msg = "Not a valid file";
    return recognizer;
  }

  recognizer->n_characters = header[1];
  recognizer->n_groups = header[2];
  recognizer->dimension = VEC_DIM_MAX;
  recognizer->downsample_threshold = header[4];

  if (recognizer->n_characters == 0 || recognizer->n_groups == 0) {
    recognizer->error_msg = "No characters in this model";
    return recognizer;
  }

  const char *cursor = recognizer->data + 5 * sizeof(unsigned int);
  recognizer->characters = (void *)cursor;

  cursor += recognizer->n_characters * sizeof(wagomu_character_info_t);
  recognizer->groups = (void *)cursor;

  recognizer->distm = (wagomu_result_t *)arena_alloc(
      &recognizer->arena, recognizer->n_characters * sizeof(wagomu_result_t),
      _Alignof(wagomu_result_t));

  unsigned int max_n_vectors = 0;
  for (unsigned int i = 0; i < recognizer->n_characters; i++)
    if (recognizer->characters[i].n_vectors > max_n_vectors)
      max_n_vectors = recognizer->characters[i].n_vectors;

  recognizer->dtw1 = (float *)arena_alloc(
      &recognizer->arena, max_n_vectors * recognizer->dimension * sizeof(float),
      _Alignof(float));
  recognizer->dtw2 = (float *)arena_alloc(
      &recognizer->arena, max_n_vectors * recognizer->dimension * sizeof(float),
      _Alignof(float));

  return recognizer;
}

// `wagomu_get_error_message` returns the last error message the
// `wagomu_recognizer_t` encountered.
const char *wagomu_get_error_message(wagomu_recognizer_t *recognizer) {
  return recognizer->error_msg;
}

// `wagomu_recognizer_start_stroke` starts a new stroke.
void wagomu_recognizer_start_stroke(wagomu_recognizer_t *recognizer) {
  if (!recognizer)
    return;

  strokes_start_stroke(&recognizer->arena, &recognizer->strokes);
}

// `wagomu_recognizer_add_stroke_point` appends a point to the latest stroke.
void wagomu_recognizer_add_stroke_point(wagomu_recognizer_t *recognizer,
                                        float x, float y) {
  if (!recognizer)
    return;

  wagomu_point_t point = {x, y};
  strokes_add_point(&recognizer->arena, &recognizer->strokes, point);
}

// `wagomu_recognize` performs a recognition on the internal strokes.
unsigned int wagomu_recognize(wagomu_recognizer_t *recognizer,
                              wagomu_result_t *results,
                              unsigned int max_results) {

  if (max_results == 0)
    return 0;

  if (!results)
    return 0;

  if (!recognizer)
    return 0;

  if (recognizer->strokes.n_strokes == 0)
    return 0;

  unsigned int total_points = recognizer->strokes.points_size;

  if (total_points < 2)
    return 0;

  // Copy strokes to processing buffer
  strokes_copy(&recognizer->arena, &recognizer->proc_strokes,
               &recognizer->strokes);

  // Pre-process strokes
  normalize_strokes(&recognizer->proc_strokes);
  downsample_strokes(&recognizer->proc_strokes,
                     recognizer->downsample_threshold);

  total_points = recognizer->proc_strokes.points_size;

  if (total_points < 2)
    return 0;

  unsigned int n_intervals = total_points - 1;

  if (recognizer->features_cap < n_intervals * recognizer->dimension) {
    recognizer->features_cap = n_intervals * recognizer->dimension;
    recognizer->features = (float *)arena_alloc(
        &recognizer->arena, recognizer->features_cap * sizeof(float),
        _Alignof(float));
  }

  // Feature extraction
  unsigned int k = 0;
  wagomu_point_t *prev;

  for (unsigned int s = 0; s < recognizer->proc_strokes.n_strokes; ++s) {
    unsigned int start = recognizer->proc_strokes.meta[s].start;
    unsigned int count = recognizer->proc_strokes.meta[s].count;

    for (unsigned int j = 0; j < count; ++j) {
      wagomu_point_t *curr = &recognizer->proc_strokes.points[start + j];

      if (prev) {
        float *c = (float *)curr;
        float *p = (float *)prev;

        for (unsigned int h = 0; h < recognizer->dimension; ++h) {
          if (h >= 2) {
            recognizer->features[k * recognizer->dimension + h] = 0.0f;
          } else {
            recognizer->features[k * recognizer->dimension + h] =
                fabsf(c[h] - p[h]);
          }
        }
        k++;
      }

      prev = curr;
    }
  }

  if (n_intervals == 0) {
    return 0;
  }

  // Matching against model
  unsigned int char_id = 0;
  unsigned int n_chars = 0;
  float *cursor;
  unsigned int n_strokes = recognizer->proc_strokes.n_strokes;

  for (unsigned int group_id = 0; group_id < recognizer->n_groups; group_id++) {
    if (recognizer->groups[group_id].n_strokes >
        (n_strokes + recognizer->window_size))
      break;

    if (recognizer->groups[group_id].n_strokes + recognizer->window_size <
        n_strokes) {
      char_id += recognizer->groups[group_id].n_chars;
      continue;
    }

    cursor = (float *)(recognizer->data + recognizer->groups[group_id].offset);
    unsigned int n_group_chars = recognizer->groups[group_id].n_chars;

    for (unsigned int i = 0; i < n_group_chars; i++) {
      recognizer->distm[n_chars].unicode =
          recognizer->characters[char_id].unicode;
      recognizer->distm[n_chars].dist =
          dtw(recognizer, recognizer->features, n_intervals, cursor,
              recognizer->characters[char_id].n_vectors);

      cursor +=
          recognizer->characters[char_id].n_vectors * recognizer->dimension;
      char_id++;
      n_chars++;
    }
  }

  qsort((void *)recognizer->distm, (unsigned int)n_chars,
        sizeof(wagomu_result_t), char_dist_cmp);

  unsigned int count = MIN(n_chars, max_results);
  if (count == 0)
    return 0;

  for (unsigned int i = 0; i < count; i++) {
    results[i] = recognizer->distm[i];
  }

  return count;
}

// `wagomu_recognizer_reset` resets the strokes back to zero.
void wagomu_recognizer_reset(wagomu_recognizer_t *recognizer) {
  if (!recognizer)
    return;

  recognizer->strokes.points_size = 0;
  recognizer->strokes.n_strokes = 0;
}

// `wagomu_recognizer_destroy` frees the arena containing the recognizer.
void wagomu_recognizer_destroy(wagomu_recognizer_t *recognizer) {
  if (!recognizer)
    return;

  wagomu_arena_t arena = recognizer->arena;
  arena_destroy(&arena);
}
