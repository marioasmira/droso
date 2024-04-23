# positive.R

#' Transform values to positive or zero
#'
#' @param x numeric value to transform
#'
#' @return transformed numeric value
#' @export
positive <- function(x) {
  y <- numeric()
  for(i in seq_along(x)){
    if(x[i] < 0)
      y <- c(y, 0)
    else
      y <- c(y, x[i])
  }
    return(y)
}
