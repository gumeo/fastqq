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

  # User is warned
  pvector <- validate_and_warn_pvector(pvector)

  # Observed and expected
  o <- rev(-log10(sort(pvector,decreasing=FALSE)))
  e <- rev(-log10(stats::ppoints(length(pvector) )))

  OEmat <- drop_dense(o, e, 1)

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

#' Internal function to prune quantiles of non-important values for
#' visualization.
#'
#' This function is not exposed, since we want to hard-code the parameters
#' for simplicity of usage.
#'
#' @param o A numeric vector of ascending sorted sample/theoretical points.
#' @param e A numeric vector of ascending sorted theoretical/sample points.
#' @param n_inp To keep track of which iteration we are in, since this is
#'              recursive.
#' @param N_hard Desired upperbound on number of points to plot.
#' @param max_iter Maximum number of rounds to try to prune points.
#' @return data.frame with o and e pruned as columns.
#' @noRd
drop_dense <- function(o, e, n_inp, N_hard = 1e4, max_iter = 10){
  if(length(o) < N_hard | n_inp > max_iter){
    return(data.frame(o=o,e=e))
  }
  mino <- min(o)
  maxo <- max(o)
  mine <- min(e)
  maxe <- max(e)
  leno <- length(o)
  lene <- length(e)

  o_width <- maxo - mino
  e_width <- maxe - mine
  distThreshold <- (o_width)/(N_hard)

  padded_o <- c(mino - o_width, o, o + o_width)
  padded_e <- c(mine - e_width, e, e + e_width)

  # Fast l1 norm
  distSurr <- (padded_o[3:(leno+2)] - padded_o[1:leno]) +
    (padded_e[3:(leno+2)] - padded_e[1:leno])

  # Find points that are surrounded by close points
  small_inds <- which(distSurr < distThreshold)

  # Only prune the even indicies
  evens <- function(x) subset(x, x %% 2 == 0)
  small_inds <- evens(small_inds)

  # Case where we cannot remove more... So prevent stack overflow.
  if(length(small_inds) == 0){
    return(drop_dense(o, e, max_iter + 1))
  }

  return(drop_dense(o[-small_inds], e[-small_inds],
                    n_inp+1))
}
