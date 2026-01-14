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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wagomu.h"

#define MAGIC_NUMBER 0x77778888
#define VEC_DIM_MAX 4
#define NORM_SIZE 1000.0f

#undef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#undef MIN3
#define MIN3(a, b, c) (MIN((a), MIN((b), (c))))

/**
 * Type definitions.
 */

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
  unsigned int size;
  unsigned int capacity;
  wagomu_point_t *points;
} wagomu_stroke_t;

struct wagomu_recognizer_s {
  const char *data;
  const wagomu_character_info_t *characters;
  const wagomu_character_group_t *groups;
  const float *strokedata;
  float *dtw1;
  float *dtw2;
  char *error_msg;
  wagomu_result_t *distm;
  wagomu_stroke_t *strokes;
  wagomu_stroke_t *proc_strokes;
  float *features;

  unsigned int n_strokes;
  unsigned int strokes_cap;
  unsigned int proc_cap;
  unsigned int features_cap;
  unsigned int n_characters;
  unsigned int n_groups;
  unsigned int dimension;
  unsigned int downsample_threshold;
  unsigned int window_size;
};

/*
 * Helper functions.
 */

// `stroke_free` frees and resets a `wagomu_stroke_t` dynamic array.
static void stroke_free(wagomu_stroke_t *stroke) {
  if (stroke->points)
    free(stroke->points);

  stroke->points = NULL;
  stroke->size = 0;
  stroke->capacity = 0;
}

// `stroke_push` pushes (and grows) a `wagomu_point_t` into a `wagomu_stroke_t`
// dynamic array.
static void stroke_push(wagomu_stroke_t *stroke, wagomu_point_t point) {
  if (stroke->size == stroke->capacity) {
    stroke->capacity = (stroke->capacity == 0) ? 8 : stroke->capacity * 2;
    stroke->points = (wagomu_point_t *)realloc(
        stroke->points, stroke->capacity * sizeof(wagomu_point_t));
  }

  stroke->points[stroke->size++] = point;
}

// `stroke_copy` deep copies (and grows) the source `wagomu_stroke_t` dynamic
// array into the `wagomu_stroke_t` dynamic array.
static void stroke_copy(wagomu_stroke_t *destination,
                        const wagomu_stroke_t *source) {
  if (source->size == 0) {
    return;
  }

  if (destination->capacity < source->size) {
    if (destination->points)
      free(destination->points);

    destination->points =
        (wagomu_point_t *)malloc(source->size * sizeof(wagomu_point_t));
    destination->capacity = source->size;
  }

  memcpy(destination->points, source->points,
         source->size * sizeof(wagomu_point_t));
  destination->size = source->size;
}

// `dist_sq` returns the euclidean distance between two `wagomu_point_t`.
static float dist_sq(wagomu_point_t point1, wagomu_point_t point2) {
  return (point1.x - point2.x) * (point1.x - point2.x) +
         (point1.y - point2.y) * (point1.y - point2.y);
}

// `normalize_strokes` resizes the strokes into a box of size
// `NORM_SIZE`x`NORM_SIZE` and centers them.
static void normalize_strokes(wagomu_stroke_t *strokes,
                              unsigned int n_strokes) {
  float min_x = FLT_MAX, max_x = -FLT_MAX;
  float min_y = FLT_MAX, max_y = -FLT_MAX;

  // Find bounding box.
  for (unsigned int i = 0; i < n_strokes; ++i) {
    wagomu_stroke_t *stroke = &strokes[i];
    for (unsigned int j = 0; j < stroke->size; ++j) {
      if (stroke->points[j].x < min_x)
        min_x = stroke->points[j].x;
      if (stroke->points[j].x > max_x)
        max_x = stroke->points[j].x;
      if (stroke->points[j].y < min_y)
        min_y = stroke->points[j].y;
      if (stroke->points[j].y > max_y)
        max_y = stroke->points[j].y;
    }
  }

  float width = max_x - min_x;
  float height = max_y - min_y;
  if (width == 0)
    width = 1;
  if (height == 0)
    height = 1;

  // Calculate scaling ratios to fit within NORM_SIZE while maintaining aspect
  // ratio.
  float ratio_w = NORM_SIZE / width;
  float ratio_h = NORM_SIZE / height;
  float ratio = (ratio_w < ratio_h) ? ratio_w : ratio_h;

  float new_width = width * ratio;
  float new_height = height * ratio;
  float dx = (NORM_SIZE - new_width) / 2.0f;
  float dy = (NORM_SIZE - new_height) / 2.0f;

  // Apply transformation
  for (unsigned int i = 0; i < n_strokes; ++i) {
    wagomu_stroke_t *stroke = &strokes[i];
    for (unsigned int j = 0; j < stroke->size; ++j) {
      stroke->points[j].x = (stroke->points[j].x - min_x) * ratio + dx;
      stroke->points[j].y = (stroke->points[j].y - min_y) * ratio + dy;
    }
  }
}

// `downsample_strokes` reduces point count based on distance to filter out
// noise or excessive detail.
static void downsample_strokes(wagomu_stroke_t *strokes, unsigned int n_strokes,
                               unsigned int downsample_threshold) {
  for (unsigned int i = 0; i < n_strokes; ++i) {
    wagomu_stroke_t *stroke = &strokes[i];
    if (stroke->size <= 1)
      continue;

    unsigned int write_idx = 0;

    write_idx++; // Always keep the first point

    wagomu_point_t last_original = stroke->points[stroke->size - 1];

    for (unsigned int j = 1; j < stroke->size; ++j) {
      wagomu_point_t last_kept = stroke->points[write_idx - 1];
      if (dist_sq(stroke->points[j], last_kept) >=
          downsample_threshold * downsample_threshold) {
        stroke->points[write_idx++] = stroke->points[j];
      }
    }

    // Ensure the last point is preserved if it's far enough from the last kept
    // point
    wagomu_point_t last_kept = stroke->points[write_idx - 1];
    if (stroke->size > 1 && dist_sq(last_original, last_kept) > 100.0f) {
      stroke->points[write_idx++] = last_original;
    }

    stroke->size = write_idx;
  }
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
wagomu_recognizer_t *wagomu_recognizer_new(const char *model_bytes,
                                           unsigned int model_size) {
  wagomu_recognizer_t *recognizer =
      (wagomu_recognizer_t *)calloc(1, sizeof(wagomu_recognizer_t));

  if (model_bytes == NULL) {
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

  recognizer->strokedata =
      (float *)(recognizer->data + (recognizer->groups)[0].offset);

  recognizer->distm = (wagomu_result_t *)malloc(recognizer->n_characters *
                                                sizeof(wagomu_result_t));

  unsigned int max_n_vectors = 0;
  for (unsigned int i = 0; i < recognizer->n_characters; i++)
    if (recognizer->characters[i].n_vectors > max_n_vectors)
      max_n_vectors = recognizer->characters[i].n_vectors;

  recognizer->dtw1 =
      (float *)malloc(max_n_vectors * recognizer->dimension * sizeof(float));
  recognizer->dtw2 =
      (float *)malloc(max_n_vectors * recognizer->dimension * sizeof(float));

  recognizer->strokes = NULL;
  recognizer->n_strokes = 0;
  recognizer->strokes_cap = 0;

  recognizer->proc_strokes = NULL;
  recognizer->proc_cap = 0;

  recognizer->features = NULL;
  recognizer->features_cap = 0;

  return recognizer;
}

// `wagomu_get_error_message` returns the last error message the
// `wagomu_recognizer_t` encountered.
const char *wagomu_get_error_message(wagomu_recognizer_t *recognizer) {
  return recognizer->error_msg;
}

// `wagomu_recognizer_start_stroke` starts a new stroke.
// It will dynamically grow the `strokes` array if it has reached
// its capacity.
void wagomu_recognizer_start_stroke(wagomu_recognizer_t *recognizer) {
  if (!recognizer)
    return;

  if (recognizer->n_strokes == recognizer->strokes_cap) {
    unsigned int old_cap = recognizer->strokes_cap;
    recognizer->strokes_cap =
        (recognizer->strokes_cap == 0) ? 4 : recognizer->strokes_cap * 2;
    recognizer->strokes = (wagomu_stroke_t *)realloc(
        recognizer->strokes, recognizer->strokes_cap * sizeof(wagomu_stroke_t));

    memset(&recognizer->strokes[old_cap], 0,
           (recognizer->strokes_cap - old_cap) * sizeof(*recognizer->strokes));
  }

  recognizer->strokes[recognizer->n_strokes].size = 0;
  recognizer->n_strokes++;
}

// `wagomu_recognizer_add_stroke_point` appends a point to the
// latest created stroke.
void wagomu_recognizer_add_stroke_point(wagomu_recognizer_t *recognizer,
                                        float x, float y) {
  if (!recognizer)
    return;

  if (recognizer->n_strokes == 0) {
    return;
  }

  wagomu_point_t point = {x, y};
  stroke_push(&recognizer->strokes[recognizer->n_strokes - 1], point);
}

// `wagomu_recognize` performs a recognition on the internal
// `strokes` dynamic array of strokes against the model loaded.
// Roughly speaking the process is a follows:
// 1. Do a deep copy of the current `strokes` dynamic array into a new dynamic
// array.
// 2. Center the strokes and normalize them to fit in a `NORM_SIZE`x`NORM_SIZE`
// box.
// 3. Downsample the amount of points in a stroke to the `downsample_threshold`
// of the model to prevent excessive detail or noise.
// 4. Run the Dynamic Time Warping algorith to get the most likely
// `max_results`.
// The heap struct returned by the function must be freed by the caller.
wagomu_prediction_t *wagomu_recognize(wagomu_recognizer_t *recognizer,
                                      unsigned int max_results) {
  if (!recognizer)
    return NULL;

  if (recognizer->n_strokes == 0)
    return NULL;

  unsigned int total_points = 0;
  for (unsigned int i = 0; i < recognizer->n_strokes; ++i)
    total_points += recognizer->strokes[i].size;

  if (total_points < 2)
    return NULL;

  // Prepare processing buffer
  if (recognizer->proc_cap < recognizer->n_strokes) {
    unsigned int old_cap = recognizer->proc_cap;
    recognizer->proc_cap = recognizer->n_strokes;
    recognizer->proc_strokes = (wagomu_stroke_t *)realloc(
        recognizer->proc_strokes,
        recognizer->proc_cap * sizeof(wagomu_stroke_t));

    memset(&recognizer->proc_strokes[old_cap], 0,
           (recognizer->proc_cap - old_cap) *
               sizeof(*recognizer->proc_strokes));
  }

  for (unsigned int i = 0; i < recognizer->n_strokes; ++i) {
    stroke_copy(&recognizer->proc_strokes[i], &recognizer->strokes[i]);
  }

  // Pre-process strokes
  normalize_strokes(recognizer->proc_strokes, recognizer->n_strokes);
  downsample_strokes(recognizer->proc_strokes, recognizer->n_strokes,
                     recognizer->downsample_threshold);

  total_points = 0;
  for (unsigned int i = 0; i < recognizer->n_strokes; ++i)
    total_points += recognizer->proc_strokes[i].size;

  if (total_points < 2)
    return NULL;

  unsigned int n_intervals = total_points - 1;

  if (recognizer->features_cap < n_intervals * recognizer->dimension) {
    recognizer->features_cap = n_intervals * recognizer->dimension;
    free(recognizer->features);
    recognizer->features =
        (float *)malloc(recognizer->features_cap * sizeof(float));
  }

  // Feature extraction
  unsigned int k = 0;
  wagomu_point_t *prev, *curr;
  int first_point = 1;

  for (unsigned int i = 0; i < recognizer->n_strokes; ++i) {
    for (unsigned int j = 0; j < recognizer->proc_strokes[i].size; ++j) {
      curr = &recognizer->proc_strokes[i].points[j];

      if (!first_point) {
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
      first_point = 0;
    }
  }

  if (n_intervals == 0) {
    return NULL;
  }

  // Matching against model
  unsigned int char_id = 0;
  unsigned int n_chars = 0;
  float *cursor;

  for (unsigned int group_id = 0; group_id < recognizer->n_groups; group_id++) {
    /* Only compare the input with templates which have
       +- window_size the same number of strokes as the input */

    if (recognizer->groups[group_id].n_strokes >
        (recognizer->n_strokes + recognizer->window_size))
      break;

    /* Use addition to avoid unsigned underflow: a < b - c  <=>  a + c < b */
    if (recognizer->groups[group_id].n_strokes + recognizer->window_size <
        recognizer->n_strokes) {
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

  /* sort the results with glibc's quicksort */
  qsort((void *)recognizer->distm, (unsigned int)n_chars,
        sizeof(wagomu_result_t), char_dist_cmp);

  unsigned int count = MIN(n_chars, max_results);

  wagomu_prediction_t *prediction = NULL;

  if (count > 0) {
    prediction = calloc(1, sizeof(wagomu_prediction_t));
    prediction->count = (unsigned int)count;
    prediction->results = malloc(count * sizeof(wagomu_result_t));
    for (unsigned int i = 0; i < count; i++) {
      prediction->results[i] = recognizer->distm[i];
    }
  }

  return prediction;
}

// `wagomu_recognizer_reset_stroke` resets the dynamic array of `strokes` back
// to zero while retaining capacity, to allow for reuse of the array without
// relocation.
void wagomu_recognizer_reset_stroke(wagomu_recognizer_t *recognizer) {
  if (!recognizer)
    return;

  for (unsigned int i = 0; i < recognizer->n_strokes; ++i) {
    recognizer->strokes[i].size = 0;
  }

  recognizer->n_strokes = 0;
}

// `wagomu_recognizer_destroy` frees the `wagomu_recognizer_t` and it's
// underlying components.
void wagomu_recognizer_destroy(wagomu_recognizer_t *recognizer) {
  if (recognizer->distm)
    free(recognizer->distm);
  if (recognizer->dtw1)
    free(recognizer->dtw1);
  if (recognizer->dtw2)
    free(recognizer->dtw2);

  if (recognizer->strokes) {
    for (unsigned int i = 0; i < recognizer->strokes_cap; ++i)
      stroke_free(&recognizer->strokes[i]);
    free(recognizer->strokes);
  }
  if (recognizer->proc_strokes) {
    for (unsigned int i = 0; i < recognizer->proc_cap; ++i)
      stroke_free(&recognizer->proc_strokes[i]);
    free(recognizer->proc_strokes);
  }
  if (recognizer->features)
    free(recognizer->features);

  free(recognizer);
}

// `wagomu_prediction_destroy` frees the prediction object and the results.
void wagomu_prediction_destroy(wagomu_prediction_t *prediction) {
  if (prediction) {
    if (prediction->results)
      free(prediction->results);
    free(prediction);
  }
}
