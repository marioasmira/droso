# find_density_peaks.R

#' Finds the density peaks in the provided vector
#'
#' @param numbers The vector for which to find the peaks.
#' @param npeaks How many peaks to output.
#' @param decreasing If the order should be decreasing or not.
#'
#' @return The density npeaks in the provided vector
#' @export
#' @importFrom stats density
find_density_peaks <- function(numbers,
                               npeaks = NULL,
                               decreasing = TRUE) {
  # Estimate the density
  d <- density(numbers)
  # Find local maxima
  local_maxima <- d$x[c(FALSE, diff(diff(d$y) >= 0) < 0)]
  # Get corresponding y values for local maxima
  y_values <- d$y[match(local_maxima, d$x)]
  # Sort y values in decreasing order
  sorted_y_values <- sort(y_values, decreasing = decreasing)
  # The second peak in terms of y values
  if (is.null(npeaks)) {
    npeaks <- length(local_maxima)
  }
  if (npeaks <= length(local_maxima)) {
    second_peak_y <- sorted_y_values[1:npeaks]
  } else {
    stop("The requested number of peaks is
    larger than the found number of peaks.")
  }
  # Corresponding x value for the second peak
  second_peak_x <-
    sort(local_maxima[which(y_values %in% second_peak_y)],
         decreasing = decreasing)
  return(second_peak_x)
}