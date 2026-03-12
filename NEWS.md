# fastqq

## fastqq 0.1.4

- Added `zero_action` parameter to `qq()`: replaces p-values of exactly zero
  with a user-supplied finite substitute (e.g. `1e-300`) and emits a warning
  stating how many values were replaced.
- New function `qqlog()`: accepts pre-computed −log₁₀(p-values) directly.
  Same signature and performance characteristics as `qq()`.
- New function `qqchisq1()`: accepts χ² test statistics (df = 1) and converts
  them to −log₁₀(p) using log-space `pchisq()` for numerical precision beyond
  `.Machine$double.xmin`.
- Added unit tests (testthat) covering the new functions.

## fastqq 0.1.3

Second cran release

- Removed bitwise or in C++ code because of warning.
- Rebuilt readme.


## fastqq 0.1.0

This is the first version released on github. The following is implemented:

- `fastqq:qq` a significantly faster alternative to `qqman::qq` for large input.
- `fastqq::qqnorm` and `fastqq::qqplot`, which are faster than their `stats` equivalents.
- `drop_dense`, which allows for retrieving the data, to use with `ggplot` or other plotting libraries.
