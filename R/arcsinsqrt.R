# arcsinsqrt.R

#' Arc sin transformation
#' @description
#' Performs the arc sin square root transformation.
#'
#' @param x Value to transform.
#' @returns The arc sin square root transformation of x.
#' @export
arcsinsqrt <- function(x) {
  return(asin(sqrt(x)))
}