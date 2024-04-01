# squarelike.R

#' Square-like function
#' @description
#' Calculates the output using a square-like shaped curve.
#'
#' @param x Point or vector to calculate the function.
#' @param l_inflection The left inflection point.
#' @param l_slope The slope of the left side of the function.
#' @param r_inflection The right inflection point.
#' @param r_slope The slope of the right side of the function.
#' @param max The maximum value for the function.
#' @returns The y value for x.
#' @export
#' @importFrom rcspline logis
squarelike <- function(
    x,
    l_inflection,
    l_slope,
    r_inflection,
    r_slope,
    max) {
  logis(x,
    x_0 = l_inflection,
    k = l_slope,
    max = max
  ) +
    logis(x,
      x_0 = r_inflection,
      k = r_slope,
      max = max
    ) -
    max
}
