# notin.R

#' Inverse of the %in% operator
#'
#' @param x A vector of values to look for in `table`.
#' @param table A vector where `x` values might be found.
#' @return A logical vector indicating if a match was not located for each
#' element of `x` in `table`.
#' @examples
#' # Check if values are in a vector
#' c(1, 2, 3) %notin% c(1, 2, 5)  # Returns: FALSE FALSE TRUE
#' @export
`%notin%` <- function(x, table) {
  return(!(x %in% table))
}