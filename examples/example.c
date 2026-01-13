#include "wagomu.h"
#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    // The model used when I was writing this example was Japanese (all) from
    // https://tegaki.github.io/.
    // If you load it, it should display U+4E00 or 一.
    fprintf(stderr, "Usage: %s <path_to_model_file>\n", argv[0]);
    return 1;
  }

  printf("%s\n", argv[1]);

  wagomu_recognizer_t *recognizer = wagomu_recognizer_new(argv[1]);

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

  wagomu_prediction_t *prediction = wagomu_recognize(recognizer, 5);

  for (unsigned int i = 0; i < prediction->count; ++i) {
    unsigned int unicode = prediction->results[i].unicode;
    float dist = prediction->results[i].dist;

    printf("  Rank %u: U+%04X (Distance: %f)\n", i + 1, unicode, dist);
  }

  wagomu_recognizer_reset_stroke(recognizer);

  wagomu_recognizer_destroy(recognizer);

  return 0;
}
