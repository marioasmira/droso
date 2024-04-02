# area_squarelike.R

#' Area under the squarelike function
#'
#' @param lower_bound Value for the lower bound of area.
#' @param upper_bound Value for the upper bound of area.
#' @param l_inflection The left inflection point.
#' @param l_slope The slope of the left side of the function.
#' @param r_inflection The right inflection point.
#' @param r_slope The slope of the right side of the function.
#' @param max The maximum value for the function.
#'
#' @return numeric value for the area under the squarelike function
#' @export
area_squarelike <- function(
    lower_bound,
    upper_bound,
    l_inflection,
    l_slope,
    r_inflection,
    r_slope,
    max
) {
  return(
    (
      antider_logis(upper_bound, max, l_inflection, l_slope) -
        antider_logis(lower_bound, max, l_inflection, l_slope)
      ) +
      (
        antider_logis(upper_bound, max, r_inflection, r_slope) -
          antider_logis(lower_bound, max, r_inflection, r_slope)
      ) -
      (max * upper_bound - max * lower_bound)
  )
}
