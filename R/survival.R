# survival.R

#' Fly survival
#' @description
#' Calculates the survival of a fly by creating a square-shaped function with
#' two logistic curves and a modified normal distribution pdf. This is a wrapper
#' for the squarelike function with the survival specific values already filled
#' in.
#'
#' @param temperature The current temperature for which to test survival
#' @param experience The phenotype (midpoint) of the fly to be tested.
#' @returns The probability of survival of the fly.
#' @export
#' @importFrom rcspline logis
survival <- function(temperature, experience) {
  squarelike(
    x = temperature,
    l_inflection = experience - 9,
    l_slope = 1.0,
    r_inflection = experience + 9,
    r_slope = -1.0,
    max = 0.8
  )
}
