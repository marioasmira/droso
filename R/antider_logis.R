# antider_logis.R

#' Anti derivative of the logistic function
#'
#' @param x Value to calculate.
#' @param max Asymptote of the logistic function.
#' @param infl Inflection point of the logistic function.
#' @param slope Slope of the logistic function.
#'
#' @return Numeric value
antider_logis <- function(
    x,
    max,
    infl,
    slope) {
  return(
    max * (x - infl) +
      (max / slope) * log(abs(1 + exp(-slope * (x - infl))))
  )
}
