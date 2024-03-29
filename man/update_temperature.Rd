% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_temperature.R
\name{update_temperature}
\alias{update_temperature}
\title{Update temperature}
\usage{
update_temperature(
  temperature,
  day,
  autocorr,
  base_env,
  amplitude,
  rand_env,
  period
)
}
\arguments{
\item{temperature}{The temperature of the previous day.}

\item{day}{The location in the sine wave to be calculated.}

\item{autocorr}{The environmental autocorrelation between days.}

\item{base_env}{The base value around which the sine wave will change.}

\item{amplitude}{The amplitude of the sine wave.}

\item{rand_env}{Amount of error to add to the environment.}

\item{period}{The period of the sine wave.}
}
\value{
Returns the value of temperature of the specified day.
}
\description{
Calculates the temperature of a specific day, taking into account the
previous
day's temperature.
}
