# survival.R

#' Fly survival
#' @description
#' Calculates the survival of a fly by creating a square-shaped function with
#' two logistic curves and a modified normal distribution pdf.
#'
#' @param temperature The current temperature for which to test survival
#' @param experience The phenotype (midpoint) of the fly to be tested.
#' @returns The probability of survival of the fly.
#' @export
#' @importFrom rcspline logis
survival <- function(temperature, experience) {
  logis(temperature,
    x_0 = experience - 9,
    k = 1.0,
    max = 0.8
  ) +
    logis(temperature,
      x_0 = experience + 9,
      k = -1.0,
      max = 0.8
    ) -
    0.8
}
