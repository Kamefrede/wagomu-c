# Wagomu C

This project is a fork of the handwritten recognition engine `wagomu` from the Tegaki project (https://github.com/tegaki/tegaki/blob/master/tegaki-engines/tegaki-wagomu/wagomu.cpp). It has been re-written in pure C to expose a smaller and simpler API.

## Motivation

The main motivation for this rewrite is to use this library as part of another project via Foreign Function Interface (FFI). Consequently, the API surface has been kept as minimal as possible while ensuring the library retains most of its original functionality.

## Differences

This implementation includes functionality that was not part of the original engine but was typically handled by external wrapper scripts. Specifically, this library handles:

- Stroke normalization
- Centering
- Downsampling

This decision was made as I personally think this responsability falls onto the library and not the caller.

## Usage

There are no precompiled libraries provided. To use this in your application, simply add `wagomu.h` and `wagomu.c` to your project source.

Please look at `examples/example.c` for an implementation example.

This library should be compatible with all existing wagomu compatible models.
