# sine_wave.R

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