# generate_reps.R

#' Digit expansion
#' @description
#' Expands a vector of tens to have all digits (19 -> 190, 191,192,...).
#'
#' @param rep_group Vector with tens to expand.
#' @returns A vector with the expanded group of numbers.
#' @export
generate_reps <- function(rep_group) {
  plot_replicate <- numeric()
  for (par in rep_group) {
    for (i in 0:9) {
      plot_replicate <- c(plot_replicate, par * 10 + i)
    }
  }
  return(plot_replicate)
}