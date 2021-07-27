
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastqq

<!-- badges: start -->
<!-- badges: end -->

`fastqq` is intended to allow users that create quantile-quantile plots
for large scale association studies, like genome wide association
studies (GWAS). In these cases, the user often plots tens to hundreds of
millions of points. Creating scatter plots with so many points is
usually not efficient, since the graphics devices store all the data,
such that the plot can be rescaled or plotted in a vector graphics
format, (where again all the data is stored).

A better and faster approach in these cases is to directly rasterize the
plot. This means that we cannot rescale the figure afterwards, and
specify the size of the plot upfront. This approach scales much better
than conventional scatter plotting with `base::plot` or
`ggplot2::ggplot` with `geom_points()`. The rasterization is achieved
with the [`scattermore` package](https://github.com/exaexa/scattermore).

See the examples below on how this compares to other styles of plotting.

**Note** that this package is inspired by the `qqman` package, which has
now [been archived](https://github.com/stephenturner/qqman). The
interface to the `qq` function should be very similar, and `fastqq::qq`
should ideally be a drop in replacement for `qqman::qq`. I created this
package, since it could take 30-60 minutes to render a single plot with
`qqman::qq`. There are other ways to inspect this type of plot faster,
but I think that `scattermore` is a good solution, since it can easily
replace existing code.

## Installation

You can install the released version of fastqq from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fastqq")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gumeo/fastqq")
```

## Example

The following is an example from very simple simulated data:

``` r
library(fastqq)
set.seed(42)
p_simulated <- runif(10000)
fastqq::qq(p_simulated)
```

<img src="man/figures/README-example-1.png" width="100%" />

We can compare the timings of creating the plots, with `qqman`.

``` r
N_test <- c(1e3,1e4,1e5,1e6)
time_method <- function(pkg_name){
  suppressPackageStartupMessages(library(pkg_name, 
                                         character.only=TRUE, quietly = TRUE))
  for(N in N_test){
    p_vec <- runif(n = N)
    print(paste0("Timing ", pkg_name, "::qq with ", 
                 N, " points"))
    tictoc::tic()
    pdf(file = NULL) # Just to prevent the plot from appearing in the readme...
    do.call("qq", list(pvector=p_vec))
    dev.off()
    tictoc::toc()  
  }
}

time_method('fastqq')
#> [1] "Timing fastqq::qq with 1000 points"
#> 0.038 sec elapsed
#> [1] "Timing fastqq::qq with 10000 points"
#> 0.02 sec elapsed
#> [1] "Timing fastqq::qq with 1e+05 points"
#> 0.048 sec elapsed
#> [1] "Timing fastqq::qq with 1e+06 points"
#> 0.411 sec elapsed
time_method('qqman')
#> [1] "Timing qqman::qq with 1000 points"
#> 0.003 sec elapsed
#> [1] "Timing qqman::qq with 10000 points"
#> 0.025 sec elapsed
#> [1] "Timing qqman::qq with 1e+05 points"
#> 0.245 sec elapsed
#> [1] "Timing qqman::qq with 1e+06 points"
#> 2.474 sec elapsed
```

So we can expect around 8-10X speedup.
