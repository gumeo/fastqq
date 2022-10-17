# fastqq

## fastqq 0.1.2

Second cran release

- Removed bitwise or in C++ code because of warning.
- Rebuilt readme.


## fastqq 0.1.0

This is the first version released on github. The following is implemented:

- `fastqq:qq` a significantly faster alternative to `qqman::qq` for large input.
- `fastqq::qqnorm` and `fastqq::qqplot`, which are faster than their `stats` equivalents.
- `drop_dense`, which allows for retrieving the data, to use with `ggplot` or other plotting libraries.
