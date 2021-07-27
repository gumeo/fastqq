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
#' @importFrom scattermore scattermoreplot
#'
#' @examples
#' \dontrun{
#' qq(stats::runif(1e6))
#' }
#' @export
qq <- function(pvector, ...) {

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
  # User has been warned...
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]

  # Observed and expected
  o = -log10(sort(pvector,decreasing=FALSE))
  e = -log10(stats::ppoints(length(pvector) ))

  def_args <- list(cex=2,
                   xlab=expression(Expected~~-log[10](italic(p))),
                   ylab=expression(Observed~~-log[10](italic(p))))
  dotargs <- list(...)
  base::tryCatch(do.call("scattermoreplot",
                   c(list(x=e, y=o),
                     def_args[!names(def_args) %in% names(dotargs)], dotargs)),
           warn=stop)

  # Add diagonal
  graphics::abline(0,1,col="red")

  # Later add shadow region.
}
