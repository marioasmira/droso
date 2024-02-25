# round_any.R

#' Round any
#' @description
#' Rounds, floors, or ceilings with provided accuracy. Copied from plyr.
#'
#' @param x Value orr vector to round.
#' @param accuracy Value with which to round.
#' @param f Function to use. Can be floor, ceiling, or round. Default is round.
#' @returns The transformed value.
# from plyr (to avoid conflicts with dplyr)
#' @export
round_any <-
  function(x, accuracy, f = round) {
    stopifnot(accuracy > 0)
    f(x / accuracy) * accuracy
  }