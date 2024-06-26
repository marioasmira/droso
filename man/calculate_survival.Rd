% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fly.R
\name{calculate_survival}
\alias{calculate_survival}
\alias{calculate_survival,fly-method}
\title{fly method to calculate survival of flies}
\usage{
calculate_survival(
  object,
  received_PE,
  temperature,
  offspring_error = NULL,
  overall_GW,
  presampled_error = NULL
)

\S4method{calculate_survival}{fly}(
  object,
  received_PE,
  temperature,
  offspring_error = NULL,
  overall_GW,
  presampled_error = NULL
)
}
\arguments{
\item{object}{An object}

\item{received_PE}{Prediction received by the mother.}

\item{temperature}{Value of temperature.}

\item{offspring_error}{Amount of error that offspring have to measure
temperature. Varies between replicates.}

\item{overall_GW}{Baseline value for GW. Varies between replicates.}

\item{presampled_error}{Pre-sampled error in cases where randomness will
be re-used}
}
\value{
A named list with the midpoint and survival according to
temperature, prediction and phenotype.
}
\description{
fly method to calculate survival of flies
}
