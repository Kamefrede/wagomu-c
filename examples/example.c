#include "../wagomu.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *read_file(const char *path, unsigned int *out_size) {
  FILE *file = fopen(path, "rb");
  if (!file)
    return NULL;

  if (fseek(file, 0, SEEK_END) != 0) {
    fclose(file);
    return NULL;
  }

  long size = ftell(file);
  if (size < 0) {
    fclose(file);
    return NULL;
  }
  rewind(file);

  char *buffer = (char *)malloc(size);
  if (!buffer) {
    fclose(file);
    return NULL;
  }

  unsigned int read = fread(buffer, 1, size, file);
  fclose(file);

  if (read != (unsigned int)size) {
    free(buffer);
    return NULL;
  }

  if (out_size)
    *out_size = size;

  return buffer;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    // The model used when I was writing this example was Japanese (all) from
    // https://tegaki.github.io/.
    // If you load it, it should display U+4E00 or 一.
    fprintf(stderr, "Usage: %s <path_to_model_file>\n", argv[0]);
    return 1;
  }

  unsigned int model_size = 0;
  const char *model = read_file(argv[1], &model_size);

  if (model == NULL) {
    fprintf(stderr, "Could not load model file %s\n", argv[1]);

    return 1;
  }

  wagomu_recognizer_t *recognizer = wagomu_recognizer_new(model, model_size);

  if (wagomu_get_error_message(recognizer) != NULL) {
    fprintf(stderr, "Error loading model: %s\n",
            wagomu_get_error_message(recognizer));

    wagomu_recognizer_destroy(recognizer);
    return 1;
  }

  wagomu_recognizer_start_stroke(recognizer);

  float x1 = 100.0f;
  float y1 = 500.f;
  float x2 = 900.f;
  float y2 = 500.f;

  float len = hypotf(x2 - x1, y2 - y1);
  int steps = (int)(len / 5.0f);
  if (steps == 0)
    steps = 1;

  for (int i = 0; i <= steps; ++i) {
    float t = (float)i / steps;
    wagomu_recognizer_add_stroke_point(recognizer, x1 + (x2 - x1) * t,
                                       y1 + (y2 - y1) * t);
  }

  unsigned int max_results = 5;
  wagomu_result_t results[5] = {0};
  unsigned int result_count =
      wagomu_recognize(recognizer, results, max_results);

  for (unsigned int i = 0; i < result_count; ++i) {
    unsigned int unicode = results[i].unicode;
    float dist = results[i].dist;

    printf("  Rank %u: U+%04X (Distance: %f)\n", i + 1, unicode, dist);
  }

  wagomu_recognizer_reset(recognizer);

  wagomu_recognizer_destroy(recognizer);

  return 0;
}
