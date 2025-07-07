#' Inverting a number around a specific value
#' @description
#' Function to invert numeric values around a specific value.
#'
#' @param x The value to be inverted.
#' @param invert_axis The axis to which x will be inverted around.
#'
#' @returns A numeric value inverted around invert_axis.
#' @export
#'
invert_value <- function(
  x,
  invert_axis
  ) {
  return(2 * invert_axis - x)
}
