# fly.R

#' Fly class
#' @description
#' Class to keep individual values together inside a single "fly".
#' Will be created with the specific values for each genomic trait and the
#' phenotype will be calculated with internal functions
#'
#' @slot replicate The simulation replicate the fly was in.
#' @slot population Population in which the fly was in.
#' @slot min_pe Minimum value for the prediction interval.
#' @slot max_pe Maximum value for the prediction interval.
#' @slot mid_pe Midpoint value for the prediction interval.
#' @slot total_pe Total length value for the prediction interval.
#' @slot precision Precision for all calculations.
#' @slot x_pe The range of used temperatures.
#' @slot pe_genes Vector of 25 values that define the PE reaction norm.
#' @slot gw_genes Vector of 5 vales that define the GW reaction norm.
#' @slot ow_genes Vector of 5 vales that define the OW reaction norm.
#' @slot pw_genes Vector of 5 vales that define the PW reaction norm.
#' @slot mean_surv Genetic value for the survival midpoint.
#' @slot PE PE reaction norm.
#' @slot GW GW reaction norm.
#' @slot OW OW reaction norm.
#' @slot PW PW reaction norm.
#' @export
setClass(
  "fly",
  slots = list(
    replicate = "numeric",
    population = "numeric",
    min_pe = "numeric",
    max_pe = "numeric",
    mid_pe = "numeric",
    total_pe = "numeric",
    precision = "numeric",
    x_pe = "numeric",
    # PE_genes should be a vector of 25 values
    pe_genes = "numeric",
    # The following genes should be a vector of 5 values each
    gw_genes = "numeric",
    ow_genes = "numeric",
    pw_genes = "numeric",
    # mean_surv is a single value
    mean_surv = "numeric",
    # Below are phenotypes
    PE = "matrix",
    GW = "numeric",
    OW = "numeric",
    PW = "numeric"
  ),
  # Default fly
  prototype = list(
    replicate = 0,
    population = 0,
    min_pe = -10,
    max_pe = 50,
    mid_pe = 20,
    total_pe = 60,
    precision = 100,
    x_pe = seq(0, 60, length.out = 100),
    # PE_genes should be a vector of 25 values
    pe_genes = c(24.4, rep(x = 0, times = 24)),
    # The following genes should be a vector of 5 values each
    pe_genes = c(3, 0, 0, 0, 0),
    ow_genes = c(-2, 0, 0, 0, 0),
    pw_genes = c(-2, 0, 0, 0, 0),
    # mean_surv is a single value
    mean_surv = 22,
    # Below are phenotypes
    PE = matrix(0, 100, 100),
    GW = rep(0, 100),
    OW = rep(0, 100),
    PW = rep(0, 100)
  )
)

#' Constructor method for fly
#'
#' @param replicate The simulation replicate the fly was in.
#' @param population Population in which the fly was in.
#' @param min_pe min_PE Minimum value for the prediction interval.
#' @param max_pe Maximum value for the prediction interval.
#' @param precision Precision for all calculations.
#' @param x_pe The range of used temperatures.
#' @param pe_genes Vector of 25 values that define the PE reaction norm.
#' @param gw_genes Vector of 5 vales that define the GW reaction norm.
#' @param ow_genes Vector of 5 vales that define the OW reaction norm.
#' @param pw_genes Vector of 5 vales that define the PW reaction norm.
#' @param mean_surv Genetic value for the survival midpoint.
#' @importFrom methods new
#' @export
fly <-
  function(replicate,
           population,
           min_pe,
           max_pe,
           precision,
           x_pe,
           pe_genes,
           gw_genes,
           ow_genes,
           pw_genes,
           mean_surv) {
    new(
      "fly",
      replicate = replicate,
      population = population,
      min_pe = min_pe,
      max_pe = max_pe,
      mid_pe = (min_pe + max_pe) / 2,
      total_pe = max_pe - min_pe,
      precision = precision,
      x_pe = seq(0, max_pe - min_pe, length.out = precision),
      pe_genes = pe_genes,
      gw_genes = gw_genes,
      ow_genes = ow_genes,
      pw_genes = pw_genes,
      mean_surv = mean_surv
    )
  }

#' @title fly method to calculate phenotype
#'
#' @param object An object
#' @param univar_matrix 1D matrix from rcspline.
#' @param bivar_matrix 2D matrix from rcspline.
#' @returns The same object but modified to have a calculated phenotype.
#' @rdname calculate_phenotype
#' @importFrom rcspline logis spline_2d spline_1d
#' @export
setGeneric("calculate_phenotype",
           function(object,
                    univar_matrix,
                    bivar_matrix) {
             standardGeneric("calculate_phenotype")
           })

#' @rdname calculate_phenotype
setMethod(
  "calculate_phenotype",
  signature = c(object = "fly"),
  definition = function(object,
                        univar_matrix,
                        bivar_matrix) {
    # PE
    object@PE <-
      logis(
        x = spline_2d(object@pe_genes, bivar_matrix),
        max = object@total_pe,
        k = 0.1,
        x_0 = object@mid_pe
      ) + object@min_pe
    # Performing the calculation beforehand instead of inside the
    # ifelse tests.
    # Assigning the values directly to the final object to make it more
    # memory efficient
    object@GW <- spline_1d(object@gw_genes, univar_matrix)
    object@OW <- spline_1d(object@ow_genes, univar_matrix)
    object@PW <- spline_1d(object@pw_genes, univar_matrix)
    # Performing the same truncation as in the original simulation
    object@GW <- 2^(object@GW)
    object@OW <- 2^(object@OW)
    object@PW <- 2^(object@PW)
    weight_sum <- object@GW + object@OW + object@PW
    object@GW <- object@GW / weight_sum
    object@OW <- object@OW / weight_sum
    object@PW <- object@PW / weight_sum

    return(object)
  }
)

#' @title fly method to retrieve prediction with two temperatures
#'
#' @param object An object
#' @param first_temp Value of temperature for first fly experience.
#' @param second_temp Value of temperature for second fly experience.
#' @returns The prediction from the reaction norm.
#' @rdname get_prediction
#' @export
setGeneric("get_prediction", function(object,
                                      first_temp,
                                      second_temp) {
  standardGeneric("get_prediction")
})

#' @rdname get_prediction
setMethod("get_prediction",
          "fly",
          function(object, first_temp, second_temp) {
            if (length(first_temp) != length(second_temp)) {
              stop("The two temperature objects don't have the same length.")
            } else {
              first_temp <- first_temp - object@min_pe
              second_temp <- second_temp - object@min_pe
              output <- numeric()
              for (i in seq_along(first_temp)) {
                output <- c(output,
                            object@PE[which.min(abs(object@x_pe - second_temp[i])),
                                      which.min(abs(object@x_pe - first_temp[i]))])
              }
              return(output)
            }
          })

#' @title fly method to calculate survival of flies
#'
#' @param object An object
#' @param received_PE Prediction received by the mother.
#' @param temperature Value of temperature.
#' @param offspring_error Amount of error that offspring have to measure
#' temperature. Varies between replicates.
#' @param overall_GW Baseline value for GW. Varies between replicates.
#' @param presampled_error Pre-sampled error in cases where randomness will
#' be re-used
#' @returns A named list with the midpoint and survival according to
#' temperature, prediction and phenotype.
#' @rdname calculate_survival
#' @importFrom stats rnorm
#' @export
setGeneric("calculate_survival", function(object,
                                          received_PE,
                                          temperature,
                                          offspring_error = NULL,
                                          overall_GW,
                                          presampled_error = NULL) {
  standardGeneric("calculate_survival")
})

#' @rdname calculate_survival
setMethod("calculate_survival",
          "fly",
          function(object,
                   received_PE,
                   temperature,
                   offspring_error = NULL,
                   overall_GW,
                   presampled_error = NULL) {
            if (is.null(presampled_error)) {
              offspring_experience <- temperature +
                rnorm(length(temperature), sd = offspring_error)
            } else if (is.null(offspring_error)) {
              offspring_experience <- temperature + presampled_error
            } else {
              stop("Both sources of error are NULL. Needs at least one.")
            }
            transformed_experience <- offspring_experience - object@min_pe
            weights_indices <- numeric()
            for (i in length(transformed_experience)) {
              weights_indices <- c(weights_indices,
                                   which.min(abs(object@x_pe - transformed_experience[i]
                                   )))
            }
            midpoint <-
              object@mean_surv * overall_GW +
              (1 - overall_GW) *
              (
                object@mean_surv * object@GW[weights_indices] +
                  offspring_experience * object@OW[weights_indices] +
                  received_PE * object@PW[weights_indices]
              )
            return(list(
              midpoint = midpoint,
              survival = survival(temperature, midpoint)
            ))
          })

#' @title fly method to retrieve information weights
#'
#' @param object An object
#' @param temperature Value of temperature.
#' @param weight Which weight to retrive. "GW", "OW", or "PW".
#' @param offspring_error Amount of error that offspring have to measure
#' temperature. Varies between replicates.
#' @param presampled_error Pre-sampled error in cases where randomness will
#' be re-used
#' @returns The information weight requested for the provided temperature.
#' @rdname get_weight
#' @importFrom stats rnorm
#' @export
setGeneric("get_weight", function(object,
                                  temperature,
                                  weight,
                                  offspring_error = NULL,
                                  presampled_error = NULL) {
  standardGeneric("get_weight")
})

#' @rdname get_weight
setMethod("get_weight",
          "fly",
          function(object,
                   temperature,
                   weight,
                   offspring_error = NULL,
                   presampled_error = NULL) {
            if (!(weight %in% c("GW", "OW", "PW"))) {
              stop("Specified weight should be \"GW\", \"OW\", or \"PW\".")
            }
            if (is.null(presampled_error)) {
              offspring_experience <- temperature +
                rnorm(length(temperature), sd = offspring_error)
            } else if (is.null(offspring_error)) {
              offspring_experience <- temperature + presampled_error
            } else {
              stop("Both sources of error are NULL. Needs at least one.")
            }
            offspring_experience <- offspring_experience - object@min_pe
            weights_indices <- numeric()
            for (i in seq_along(offspring_experience)) {
              weights_indices <- c(weights_indices,
                                   which.min(abs(object@x_pe - offspring_experience[i]
                                   )))
            }
            if (weight == "GW") {
              return(object@GW[weights_indices])
            } else if (weight == "OW") {
              return(object@OW[weights_indices])
            } else if (weight == "PW") {
              return(object@PW[weights_indices])
            }
          })
