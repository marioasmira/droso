#' Update temperature
#' @description
#' Calculates the temperature of a specific day, taking into account the
#' previous
#' day's temperature.
#'
#' @param temperature The temperature of the previous day.
#' @param day The location in the sine wave to be calculated.
#' @param autocorr The environmental autocorrelation between days.
#' @param base_env The base value around which the sine wave will change.
#' @param amplitude The amplitude of the sine wave.
#' @param rand_env Amount of error to add to the environment.
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
           rand_env,
           period) {
    stopifnot(autocorr >= 0 && autocorr <= 1)
    previous_temperature <- temperature
    return((
      autocorr *
        previous_temperature +
        (1 - autocorr) *
          sine_wave(base_env, amplitude, day, period)
    ) + rnorm(n = 1, mean = 0, sd = rand_env))
  }
