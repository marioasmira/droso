# custom_cor.R

#' Bootstrap correlation
#' @description
#' Calculate correlations from a dataset with two columns. Used for bootstrap.
#'
#' @param d The dataset object to use.
#' @param i The row in the dataset to calculate the correlation.
#' @returns The correlation for the selected row.
#' @export
#' @importFrom stats cor
custom_cor <- function(d, i) {
  return(cor(d[i, 1], d[i, 2], use = "complete.obs"))
}