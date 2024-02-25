# skew_norm_pdf.R

#' Normal distribution pdf
#' @description
#' Calculate the normal distribution probability density function at a specific
#' x value.
#'
#' @param x Numeric value for which to calculate the normal pdf.
#' @returns The normal pdf for x.
normal_pdf <- function(x) {
    exp(-(x^2) / 2) / sqrt(2 * pi)
}

#' Normal distribution cdf
#' @description
#' Calculate the normal distribution cumulative density function at a specific
#' x value.
#'
#' @param x Numeric value for which to calculate the normal cdf.
#' @returns The normal cdf for x.
#' @importFrom statip erf
normal_cdf <- function(x) {
    return((1 + erf(x / sqrt(2.0))) * 0.5)
}

#' Skew normal distribution pdf
#' @description
#' Calculate the skew normal distribution probability density function at a
#' specific x value.
#'
#' @param x The x value.
#' @param mean The mean of the normal distribution.
#' @param sd The standard deviation of the normal distribution.
#' @param skew The skew parameter for a skewed normal distribution.
#' @returns The normal pdf for x.
#' @export
skew_normal_pdf <- function(x, mean, sd, skew) {
    stopifnot(sd > 0)
    value <- (x - mean) / sd
    pdf <- normal_pdf(value)
    cdf <- normal_cdf(value * skew)
    return((2 / sd) * pdf * cdf)
}
