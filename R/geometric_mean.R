# geometric_mean.R

#' Geometric mean
#' @description
#' Calculates the geometric mean of a vector of values.
#'
#' @param x The vector for which to calculate the geometric mean.
#' @param na_rm If NAs should be removed from the calculation. Default is TRUE.
#' @returns The geometric mean of the vector.
#' @export
geometric_mean <- function(x, na_rm = TRUE) {
    stopifnot(all(x > 0))
  exp(mean(log(x), na.rm = na_rm))
}