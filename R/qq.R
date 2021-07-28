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
#' @importFrom scattermore scattermoreplot
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
  # qqman old example:
  # https://github.com/stephenturner/qqman/blob/v0.0.0/qqman.r#L335
}

drop_dense <- function(o,e,n_inp){
  # Prune down to under 10K
  N_hard <- 1e4
  if(length(o) < N_hard | n_inp > 10){
    print(paste0('Pruned down to ',length(o),' points for qq plot'))
    return(cbind(o,e))
  }
  mino <- min(o)
  maxo <- max(o)
  mine <- min(e)
  maxe <- max(e)

  distThreshold <- (maxo - mino)/(N_hard)

  distLo <- abs(o - c(mino - 100, o[1:(length(o)-1)]))
  distRo <- abs(o - c(o[2:length(o)], maxo + 100))
  distLe <- abs(e - c(mine - 100, e[1:(length(e)-1)]))
  distRe <- abs(e - c(e[2:length(e)], maxe + 100))

  distSurr <- distLo+distRo+distLe+distRe
  small_inds <- which(distSurr < distThreshold)

  evens <- function(x) subset(x, x %% 2 == 0)
  small_inds <- evens(small_inds)
  if(length(small_inds) == 0){
    return(drop_dense(o, e,11))
  }

  return(drop_dense(o[-small_inds], e[-small_inds],n_inp+1))
}

#' Creates a Q-Q plot, with redundant points removed
#'
#' Creates a quantile-quantile plot from p-values from a GWAS study. We compare
#' the data quantile with a theoretical quantile from a uniform distribution.
#' This code is mostly adapted from the \code{qqman} package, but improved
#' for speed. Here we drop redundant points, where they are so dense in the plot
#' that we would not see them anyways. This makes the plot much faster.
#'
#' @param pvector A numeric vector of p-values.
#' @param ... Other arguments passed to \code{plot()}
#'
#' @keywords visualization qq qqplot
#'
#' @examples
#' \dontrun{
#' qq_drop_dense(stats::runif(1e6))
#' }
#' @export
qq_drop_dense <- function(pvector, ...) {

  # User is warned
  pvector <- validate_and_warn_pvector(pvector)

  # Observed and expected
  o = -log10(sort(pvector,decreasing=TRUE))
  e = rev(-log10(stats::ppoints(length(pvector) )))

  OEmat <- drop_dense(o,e,1)

  o <- OEmat[,1]
  e <- OEmat[,2]

  # Use base plotting
  def_args <- list(pch=20, xlim=c(0, max(e)), ylim=c(0, max(o)),
                   xlab=expression(Expected~~-log[10](italic(p))),
                   ylab=expression(Observed~~-log[10](italic(p)))
  )
  dotargs <- list(...)

  tryCatch(do.call("plot", c(list(x=e, y=o), def_args[!names(def_args) %in% names(dotargs)], dotargs)), warn=stop)

  # Add diagonal
  abline(0,1,col="red")

  # Later add shadow region.
}
