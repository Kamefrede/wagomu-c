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
 */

#ifndef WAGOMU_H
#define WAGOMU_H

#ifdef __cplusplus
extern "C" {
#endif

// `wagomu_result_t` holds the a prediction result with unicode codepoint and
// distance from the provided strokelist.
typedef struct {
  unsigned int unicode;
  float dist;
} wagomu_result_t;

// `wagomu_prediction_t` holds an array of predictions.
// It is the responsibility of the one who holds this struct to free it.
typedef struct {
  unsigned int count;
  wagomu_result_t *results;
} wagomu_prediction_t;

/**
 * A recognizer tracks a list of strokes, which in turn are a dynamic array of
 * points. Whenever `wagomu_recognizer_start_stroke` is called, a new stroke is
 * appended and all subsequent calls to `wagomu_recognizer_add_stroke_point`
 * will add a point to the latest stroke.
 */

// `wagomu_recognizer_s` is an opaque handle to the internal
// `wagomu_recognizer_t`
typedef struct wagomu_recognizer_s wagomu_recognizer_t;

// `wagomu_recognizer_new` creates a new recognizer and loads the
// model.
// It will always return a non-null pointer, even in case of an error in
// initialization.
// To check if an error occured use `wagomu_get_error_message`.
wagomu_recognizer_t *wagomu_recognizer_new(const char *model_bytes,
                                           unsigned int model_size);
// `wagomu_get_error_message` returns the latest error that the recognizer
// encountered or `NULL` in case of no error.
const char *wagomu_get_error_message(wagomu_recognizer_t *recognizer);

// `wagomu_recognizer_start_stroke` starts a new stroke.
// All subsequent calls of `wagomu_recognizer_add_stroke_point` will append the
// points to the latest created stroke.
void wagomu_recognizer_start_stroke(wagomu_recognizer_t *recognizer);
// `wagomu_recognizer_add_stroke_point` adds a point to the current stroke.
void wagomu_recognizer_add_stroke_point(wagomu_recognizer_t *recognizer,
                                        float x, float y);
// `wagomu_recognize` performs recognition on the accumulated strokes.
// It returns a struct with up to `max_results` predictions.
// It is the responsibility of the caller of this function to free the
// prediction.
wagomu_prediction_t *wagomu_recognize(wagomu_recognizer_t *recognizer,
                                      unsigned int max_results);
// `wagomu_recognizer_reset_stroke` resets the current array of strokes
// allowing for a new character to be recognized.
void wagomu_recognizer_reset_stroke(wagomu_recognizer_t *recognizer);
// `wagomu_prediction_destroy` frees the prediction struct and the array of
// results.
void wagomu_prediction_destroy(wagomu_prediction_t *prediction);
// `wagomu_recognizer_destroy` frees the recognizer struct and its underlying
// elements.
void wagomu_recognizer_destroy(wagomu_recognizer_t *recognizer);

#ifdef __cplusplus
}
#endif

#endif
