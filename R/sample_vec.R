# sample_vec.R

#' Returns a random sample from a vector
#'
#' @param x Vector to sample from
#' @param ... Arguments for the sample function
#'
#' @return Returns sampled values from the supplied vector
#' @export
sample_vec <- function(x, ...) {
  return(x[sample(length(x), ...)])
}