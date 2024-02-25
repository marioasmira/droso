# droso_utils.R

#' Sine wave
#' @description
#' Calculates the sine value at point x. Incorporates a base value, amplitude,
#' and period.
#'
#' @param base_env The base value around which the sine wave will change.
#' @param amplitude The amplitude of the sine wave.
#' @param day The location in the sine wave to be calculated.
#' @param period The period of the sine wave.
#' @returns Returns the value of the sine wave on the specified day.
#' @export
sine_wave <- function(base_env, amplitude, day, period) {
  return(base_env + amplitude * sin(day * period))
}

#' Update temperature
#' @description
#' Calculates the temperature of a specific day, taking into account the previous
#' day's temperature.
#'
#' @param temperature The temperature of the previous day.
#' @param day The location in the sine wave to be calculated.
#' @param autocorr The environmental autocorrelation between days.
#' @param base_env The base value around which the sine wave will change.
#' @param amplitude The amplitude of the sine wave.
#' @param period The period of the sine wave.
#' @returns Returns the value of temperature of the specified day.
#' @export
#' @importFrom stats rnorm
update_temperature <-
  function(temperature,
           day,
           autocorr,
           base_env,
           amplitude,
           period) {
    previous_temperature = temperature
    
    return(
      autocorr * previous_temperature + (1 - autocorr) * sine_wave(base_env, amplitude, day, period)
    ) + rnorm(1)
  }

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
                  max = 0.8) +
    logis(temperature,
                    x_0 = experience + 9,
                    k = -1.0,
                    max = 0.8) -
    0.8 +
    exp(-((temperature - experience) ^ 2) / (2.0 * (9 * 2) ^ 2)) / (9 * 2 * sqrt(2 *
                                                                                   pi))
}

#' Geometric mean
#' @description
#' Calculates the geometric mean of a vector of values.
#'
#' @param x The vector for which to calculate the geometric mean.
#' @param na.rm If NAs should be removed from the calculation. Default is TRUE.
#' @returns The geometric mean of the vector.
#' @export
geometric_mean <- function(x, na.rm = T) {
  exp(mean(log(x), na.rm = na.rm))
}

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
    f(x / accuracy) * accuracy
  }

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

#' Standard error of the mean
#' @description
#' Calculates the standard error of the mean.
#'
#' @param x A vector.
#' @returns The standard error of the mean from the vector.
#' @export
#' @importFrom stats sd
std_mean <- function(x) {
  sd(x) / sqrt(length(x))
}

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

#' Arc sin transformation
#' @description
#' Performs the arc sin square root transformation.
#'
#' @param x Value to transform.
#' @returns The arc sin square root transformation of x.
#' @export
arcsinsqrt <- function(x) {
  return(asin(sqrt(x)))
}

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

#' Combine log files into a single file and delete old files
#'
#' @param log_directory directory where to find the log files
#'
#' @return Nothing
#' @export
combine_logs <- function(log_directory) {
  # Get the names of the log files
  log_files <- list.files(path = log_directory, pattern = "log_")
  # Combine the log files
  log <- unlist(lapply(paste0(log_directory, log_files), readLines))
  # Write the combined log to a file
  writeLines(log,
             con = paste0(
               log_directory,
               as.character(Sys.Date()),
               "_combined_log.txt"
             ))

  # Delete the separate log files
  file.remove(paste0(log_directory, log_files))
}