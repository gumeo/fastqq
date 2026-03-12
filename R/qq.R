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
#' @param zero_action A numeric value to substitute for p-values of exactly
#'   zero before plotting. If \code{NULL} (default), zero p-values are treated
#'   like any other invalid value and silently excluded. If a numeric value is
#'   provided (e.g. \code{1e-300}), all zero p-values are replaced with that
#'   value and a warning is emitted stating how many were replaced.
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
#'
#' # Handle p-values of zero by substituting a small finite value
#' pvec <- c(stats::runif(1e4), 0, 0)
#' qq(pvec, zero_action = 1e-300)
#' @export
qq <- function(pvector, zero_action = NULL, ...) {
  # Hard coded limit, for when we switch to pruning
  N_hard <- 10000

  # Replace zero p-values if requested
  if (!is.null(zero_action)) {
    n_zeros <- sum(pvector == 0, na.rm = TRUE)
    if (n_zeros > 0) {
      warning(paste0(n_zeros, " p-value(s) equal to zero replaced with ",
                     zero_action, "."))
      pvector[pvector == 0] <- zero_action
    }
  }

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
}

#' Creates a Q-Q plot from pre-computed -log10(p-values)
#'
#' Accepts \eqn{-\log_{10}(p)} values directly as input. This is useful when
#' the caller has already transformed their p-values, or when higher numerical
#' precision is required before passing values to the plotting layer. Produces
#' the same style of plot and uses the same fast pruning algorithm as
#' \code{\link{qq}}.
#'
#' @param log10_pvector A numeric vector of \eqn{-\log_{10}(p)} values. Values
#'   must be positive and finite; non-positive, infinite, NA, and NaN values
#'   are excluded (consistent with \code{qq} excluding p-values \eqn{\geq 1}
#'   or \eqn{\leq 0}).
#' @param ... Other arguments passed to \code{plot()}
#'
#' @keywords visualization qq qqplot
#'
#' @import graphics
#'
#' @return No return value, called for plotting side effects.
#'
#' @examples
#' pvec <- stats::runif(1e5)
#' qqlog(-log10(pvec))
#' @export
qqlog <- function(log10_pvector, ...) {
  N_hard <- 10000
  if (length(log10_pvector) <= N_hard) {
    log10_pvector <- log10_pvector[
      !is.na(log10_pvector) & !is.nan(log10_pvector) &
      is.finite(log10_pvector) & log10_pvector > 0
    ]
    o <- sort(log10_pvector, decreasing = TRUE)
    e <- -log10(stats::ppoints(length(o)))
  } else {
    OEmat <- drop_dense_qqlog(log10_pvector, N_hard)
    o <- OEmat[, 1]
    e <- OEmat[, 2]
  }

  def_args <- list(pch = 20, xlim = c(0, max(e)), ylim = c(0, max(o)),
                   xlab = expression(Expected~~-log[10](italic(p))),
                   ylab = expression(Observed~~-log[10](italic(p))))
  dotargs <- list(...)
  base::tryCatch(do.call("plot",
                   c(list(x = e, y = o),
                     def_args[!names(def_args) %in% names(dotargs)], dotargs)),
           warn = stop)

  abline(0, 1, col = "red")
}

#' Creates a Q-Q plot from chi-squared statistics
#'
#' Accepts \eqn{\chi^2} test statistics as input and converts them to
#' \eqn{-\log_{10}(p)} values assuming 1 degree of freedom. The conversion
#' uses log-space computation via \code{pchisq(..., log.p = TRUE)}, which
#' avoids floating-point underflow for very large test statistics and is
#' numerically precise well beyond \code{.Machine$double.xmin}. Produces the
#' same style of plot as \code{\link{qq}} and \code{\link{qqlog}}.
#'
#' @param chisq_vector A numeric vector of \eqn{\chi^2} statistics
#'   (non-negative). Negative, infinite, NA, and NaN values are excluded.
#' @param ... Other arguments passed to \code{plot()}
#'
#' @keywords visualization qq qqplot
#'
#' @import graphics
#'
#' @return No return value, called for plotting side effects.
#'
#' @examples
#' chisq_vals <- stats::rchisq(1e5, df = 1)
#' qqchisq1(chisq_vals)
#' @export
qqchisq1 <- function(chisq_vector, ...) {
  N_hard <- 10000
  if (length(chisq_vector) <= N_hard) {
    chisq_vector <- chisq_vector[
      !is.na(chisq_vector) & !is.nan(chisq_vector) &
      is.finite(chisq_vector) & chisq_vector >= 0
    ]
    # Log-space conversion avoids underflow for extreme chi-sq values
    log10_pvals <- -stats::pchisq(chisq_vector, df = 1,
                                  lower.tail = FALSE, log.p = TRUE) / log(10)
    log10_pvals <- log10_pvals[log10_pvals > 0]  # exclude chi-sq = 0 (p = 1)
    o <- sort(log10_pvals, decreasing = TRUE)
    e <- -log10(stats::ppoints(length(o)))
  } else {
    OEmat <- drop_dense_chisq1(chisq_vector, N_hard)
    o <- OEmat[, 1]
    e <- OEmat[, 2]
  }

  def_args <- list(pch = 20, xlim = c(0, max(e)), ylim = c(0, max(o)),
                   xlab = expression(Expected~~-log[10](italic(p))),
                   ylab = expression(Observed~~-log[10](italic(p))))
  dotargs <- list(...)
  base::tryCatch(do.call("plot",
                   c(list(x = e, y = o),
                     def_args[!names(def_args) %in% names(dotargs)], dotargs)),
           warn = stop)

  abline(0, 1, col = "red")
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
