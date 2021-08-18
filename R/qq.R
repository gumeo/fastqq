#  File R/qq.R
#  Part of the R package fastqq.
#
#  Copyright (C) 2021 Gudmundur Einarsson
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

#' Clean values in p-value array
#'
#' Helper function to sanitize and print warnings about unwanted values. Helper
#' function, which is not exposed.
#'
#' @param pvector A numeric vector of p-values.
#' @return Vector where invalid values have been removed.
#' @noRd
clean_pvector <- function(pvector){
  # Remove
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) &
                       !is.null(pvector) & is.finite(pvector) &
                       pvector<1 & pvector>0]
  return(pvector)
}

#' Creates a Q-Q plot
#'
#' Creates a quantile-quantile plot from p-values from an association study,
#' e.g. a genome wide association study (GWAS). We compare
#' the data quantile with a theoretical quantile from a uniform distribution.
#' This code is mostly adapted from the \code{qqman} package, but improved
#' for speed. A graph with a hundred million points should only take a few
#' seconds to generate.
#'
#' @param pvector A numeric vector of p-values.
#' @param ... Other arguments passed to \code{plot()}
#'
#' @keywords visualization qq qqplot
#'
#' @import graphics
#'
#' @return No return value, called for plotting side effects.
#'
#' @examples
#' qq(stats::runif(1e6))
#' @export
qq <- function(pvector, ...) {
  # Hard coded limit, for when we switch to pruning
  N_hard <- 10000
  if(length(pvector) <= N_hard){
    pvector <- clean_pvector(pvector)
    # Observed and expected
    o <- -log10(sort(pvector,decreasing=FALSE))
    e <- -log10(stats::ppoints(length(pvector)))
  }else{
    # All in cpp
    OEmat <- drop_dense_qq(pvector, N_hard)
    o <- OEmat[,1]
    e <- OEmat[,2]
  }

  def_args <- list(pch=20, xlim=c(0, max(e)), ylim=c(0, max(o)),
                   xlab=expression(Expected~~-log[10](italic(p))),
                   ylab=expression(Observed~~-log[10](italic(p)))
  )
  dotargs <- list(...)
  base::tryCatch(do.call("plot",
                   c(list(x=e, y=o),
                     def_args[!names(def_args) %in% names(dotargs)], dotargs)),
           warn=stop)

  # Add diagonal
  abline(0,1,col="red")

  # Later add shadow region.
  # qqman old example:
  # https://github.com/stephenturner/qqman/blob/v0.0.0/qqman.r#L335
}

#' Creates a Q-Q plot
#'
#' Faster alternative to \code{stats::qqplot()}. For more than 1e5 points
#' we remove excess points, that would not be visible in the plot, since the
#' points are so close.
#'
#' @param x First sample for \code{qqplot}.
#' @param y Second sample for \code{qqplot}.
#' @param xlab x label for plot.
#' @param ylab y label for plot.
#' @param plot.it Should the plot be created.
#' @param ... Other arguments passed to \code{plot()}
#'
#' @return list with sorted samples, interpolated to be same size.
#'
#' @keywords visualization qq qqplot
#'
#' @import graphics
#'
#' @examples
#' qqplot(stats::runif(1e6),stats::runif(1e6))
#' @export
qqplot <- function(x, y, plot.it = TRUE,
                   xlab = deparse1(substitute(x)),
                   ylab = deparse1(substitute(y)), ...) {
  sx_orig <- base::sort(x,decreasing = TRUE)
  sy_orig <- base::sort(y,decreasing = TRUE)
  lenx <- base::length(sx_orig)
  leny <- base::length(sy_orig)
  if (leny < lenx)
    sx_orig <- stats::approx(1L:lenx, sx_orig, n = leny)$y
  if (leny > lenx)
    sy_orig <- stats::approx(1L:leny, sy_orig, n = lenx)$y

  # Only difference from stats::qqplot
  XYmat <- drop_dense(sx_orig,sy_orig)
  # Make it faster by dropping excess points for plotting
  sx <- XYmat[,1]
  sy <- XYmat[,2]
  if (plot.it)
    plot(sx, sy, xlab = xlab, ylab = ylab, ...)
  # Return all values to be consistent with prior usage
  invisible(list(x = sx_orig, y = sy_orig))
}

#' Internal function to prune quantiles of non-important values for
#' visualization.
#'
#' This function is not exposed, since we want to hard-code the parameters
#' for simplicity of usage.
#'
#' @param x A numeric vector of sample/theoretical points.
#' @param y A numeric vector of theoretical/sample points.
#' @param N_hard Desired upper bound on the number of points to plot.
#' @return data.frame with o and e pruned as columns.
#' @export
drop_dense <- function(x, y, N_hard = 1e4){
  x <- base::sort(x, decreasing = TRUE)
  y <- base::sort(y, decreasing = TRUE)
  if(base::length(x) < N_hard){
    return(data.frame(x=x,y=y))
  }
  return(drop_dense_internal(x,y,N_hard))
}

#' Creates a Q-Q plot for comparing with normal quantiles
#'
#' Faster alternative to \code{stats::qqnorm()}. For more than 1e5 points
#' we remove excess points, that would not be visible in the plot, since the
#' points are so close. Otherwise this should work exactly the same, and the
#' code is mostly adapted from \code{stats::qqnorm()}. This code produces
#' more lightweight plots for excessive amounts of data.
#'
#' @param y sample, to compare to normal quantiles.
#' @param ylim graphical limits.
#' @param main Plot title.
#' @param xlab X label.
#' @param ylab Y label.
#' @param datax logical. Should data values be on x-axis?
#' @param plot.it Should the plot be created.
#' @param ... Other arguments passed to \code{plot()}
#'
#' @return \code{data.frame} with sorted sample and normal quantiles, \code{NA}
#' values are excluded.
#'
#' @keywords visualization qq qqplot
#'
#' @import graphics
#'
#' @examples
#' qqnorm(stats::rnorm(1e6))
#' @export
qqnorm <- function(y, ylim, main = "Normal Q-Q Plot",
                   xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                   plot.it = TRUE, datax = FALSE, ...)
{
  # Let's remove NA, and that is the default behavior here.
  # Sorting removes NA
  y <- base::sort(y, decreasing = TRUE)
  if(0 == (n <- base::length(y)))
    stop("y is empty or has only NAs")
  if (plot.it && missing(ylim))
    ylim <- range(y)
  x <- rev(stats::qnorm(stats::ppoints(n)))
  # Only difference from stats::qqplot
  XYmat <- drop_dense(x,y)
  # Make it faster by dropping excess points for plotting
  sx <- XYmat[,1]
  sy <- XYmat[,2]
  if(plot.it){
    if (datax){
      plot(sy, sx, main = main, xlab = ylab, ylab = xlab, xlim = ylim, ...)
    }else{
      plot(sx, sy, main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
    }
  }
  # Order is now inconsistent with stats::qqnorm
  invisible(if(datax) list(x = y, y = x) else list(x = x, y = y))
}
