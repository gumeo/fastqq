#' Validate values in p-value array
#'
#' Helper function to sanitize and print warnings about unwanted values. Helper
#' function, which is not exposed.
#'
#' @param pvector A numeric vector of p-values.
#' @return Vector where invalid values have been removed.
#' @noRd
validate_and_warn_pvector <- function(pvector){
  # Check for sensible input
  if (!is.numeric(pvector)) stop("Input must be numeric.")
  if(base::any(base::is.na(pvector))){
    warning("Input to qq contains one or more NA values, these will be discarded.")
  }
  if(base::any(base::is.nan(pvector))){
    warning("Input to qq contains one or more NaN values, these will be discarded.")
  }
  if(base::any(base::is.null(pvector))){
    warning("Input to qq contains one or more NULL values, these will be discarded.")
  }
  if(base::any(!base::is.finite(pvector))){
    warning("Input to qq contains one or more non finite values, these will be discarded.")
  }
  if(base::any(pvector < 0.0)){
    warning("Input to qq contains one or more negative values, these will be discarded.")
  }
  if(base::any(pvector > 1.0)){
    warning("Input to qq contains one or more values larger than one, these will be discarded.")
  }
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
  return(pvector)
}

#' Creates a Q-Q plot
#'
#' Creates a quantile-quantile plot from p-values from a GWAS study. We compare
#' the data quantile with a theoretical quantile from a uniform distribution.
#' This code is mostly adapted from the \code{qqman} package, but improved
#' for speed. A graph with a hundred million points should only take a few
#' seconds to generate. Note that the graph is rasterised, so you should
#' specify the size of the plot in pixels, (that you plan to export), in the
#' input.
#'
#' @param pvector A numeric vector of p-values.
#' @param ... Other arguments passed to \code{plot()}
#'
#' @keywords visualization qq qqplot
#'
#' @import graphics
#'
#' @examples
#' \dontrun{
#' qq(stats::runif(1e6))
#' }
#' @export
qq <- function(pvector, ...) {

  # User is warned, takes 5 seconds for 1e8 points...
  #pvector <- validate_and_warn_pvector(pvector)

  # Observed and expected
  #o <- -log10(sort(pvector,decreasing=FALSE))
  #e <- -log10(stats::ppoints(length(pvector)))

  OEmat <- drop_dense_qq(pvector)

  o <- OEmat[,1]
  e <- OEmat[,2]

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
#' @keywords visualization qq qqplot
#'
#' @import graphics
#'
#' @examples
#' \dontrun{
#' qqplot(stats::runif(1e6),stats::runif(1e6))
#' }
#' @export
qqplot <- function(x, y, plot.it = TRUE,
                   xlab = deparse1(substitute(x)),
                   ylab = deparse1(substitute(y)), ...) {
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx)
    sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx)
    sy <- approx(1L:leny, sy, n = lenx)$y

  # Only difference from stats::qqplot
  if(length(sx) > 1e5){
    XYmat <- drop_dense(sx,sy)
    # Make it faster by dropping excess points for plotting
    sx <- XYmat[,1]
    sy <- XYmat[,2]
  }
  if (plot.it)
    plot(sx, sy, xlab = xlab, ylab = ylab, ...)
  invisible(list(x = sx, y = sy))
}

#' Internal function to prune quantiles of non-important values for
#' visualization.
#'
#' This function is not exposed, since we want to hard-code the parameters
#' for simplicity of usage.
#'
#' @param o A numeric vector of descending sorted sample/theoretical points.
#' @param e A numeric vector of descending sorted theoretical/sample points.
#' @param N_hard Desired upperbound on number of points to plot.
#' @return data.frame with o and e pruned as columns.
#' @noRd
drop_dense <- function(o, e, N_hard = 1e4){
  if(length(o) < N_hard){
    return(data.frame(o=o,e=e))
  }
  return(drop_dense_internal(o,e,N_hard))
}
