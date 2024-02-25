# std_mean.R

#' Standard error of the mean
#' @description
#' Calculates the standard error of the mean.
#'
#' @param x A vector.
#' @returns The standard error of the mean from the vector.
#' @export
#' @importFrom stats sd
std_mean <- function(x) {
  sd(x) / sqrt(length(x))
}