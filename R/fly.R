# fly.R

#' Fly class
#' @description
#' Class to keep individual values together inside a single "fly".
#' Will be created with the specific values for each genomic trait and the
#' phenotype will be calculated with internal functions
#'
#' @slot replicate The simulation replicate the fly was in.
#' @slot population Population in which the fly was in.
#' @slot max_PE Maximum value for the prediction interval.
#' @slot precision Precision for all calculations.
#' @slot x_PE The range of used temperatures.
#' @slot PE_genes Vector of 25 values that define the PE reaction norm.
#' @slot GW_genes Vector of 5 vales that define the GW reaction norm.
#' @slot OW_genes Vector of 5 vales that define the OW reaction norm.
#' @slot PW_genes Vector of 5 vales that define the PW reaction norm.
#' @slot mean_surv Genetic value for the survival midpoint.
#' @slot fecundity_genes Vector of 3 values that define the mean,
#' sd and skew for fecundity.
#' @slot PE PE reaction norm.
#' @slot GW GW reaction norm.
#' @slot OW OW reaction norm.
#' @slot PW PW reaction norm.
#' @slot egg_laying Reaction for egg laying.
#' @export
setClass(
  "fly",
  slots = list(
    replicate = "numeric",
    population = "numeric",
    min_PE = "numeric",
    max_PE = "numeric",
    precision = "numeric",
    x_PE = "numeric",
    # PE_genes should be a vector of 25 values
    PE_genes = "numeric",
    # The following genes should be a vector of 5 values each
    GW_genes = "numeric",
    OW_genes = "numeric",
    PW_genes = "numeric",
    # mean_surv is a single value
    mean_surv = "numeric",
    # fecundity_genes should be a vector with 3 values:
    # mean_eggs, sd_eggs, and skew
    fecundity_genes = "numeric",
    # Below are phenotypes
    PE = "matrix",
    GW = "numeric",
    OW = "numeric",
    PW = "numeric",
    egg_laying = "numeric"
  ),
  # Default fly
  prototype = list(
    replicate = 0,
    population = 0,
    min_PE = -20,
    max_PE = 30,
    precision = 100,
    x_PE = seq(min_PE, max_PE, length.out = precision),
    # PE_genes should be a vector of 25 values
    PE_genes = rep(x = 0, times = 25),
    # The following genes should be a vector of 5 values each
    GW_genes = c(1, 0, 0, 0, 0),
    OW_genes = c(1, 0, 0, 0, 0),
    PW_genes = c(1, 0, 0, 0, 0),
    # mean_surv is a single value
    mean_surv = 25,
    # fecundity_genes should be a vector with 3 values:
    # mean_eggs, sd_eggs, and skew
    fecundity_genes = c(25, 15, -5),
    # Below are phenotypes
    PE = matrix(0, 100, 100),
    GW = rep(0, 100),
    OW = rep(0, 100),
    PW = rep(0, 100),
    egg_laying = rep(0, 100)
  )
)

#' Constructor method for fly
#'
#' @param replicate The simulation replicate the fly was in.
#' @param population Population in which the fly was in.
#' @param max_PE Maximum value for the prediction interval.
#' @param precision Precision for all calculations.
#' @param x_PE The range of used temperatures.
#' @param PE_genes Vector of 25 values that define the PE reaction norm.
#' @param GW_genes Vector of 5 vales that define the GW reaction norm.
#' @param OW_genes Vector of 5 vales that define the OW reaction norm.
#' @param PW_genes Vector of 5 vales that define the PW reaction norm.
#' @param mean_surv Genetic value for the survival midpoint.
#' @param fecundity_genes Vector of 3 values that define the mean,
#' sd and skew for fecundity.
#' @importFrom methods new
#' @export
fly <-
  function(replicate,
           population,
           min_PE,
           max_PE,
           precision,
           x_PE,
           PE_genes,
           GW_genes,
           OW_genes,
           PW_genes,
           mean_surv,
           fecundity_genes) {
    new(
      "fly",
      replicate = replicate,
      population = population,
      min_PE = min_PE,
      max_PE = max_PE,
      precision = precision,
      x_PE = seq(0, max_PE, length.out = precision),
      PE_genes = PE_genes,
      GW_genes = GW_genes,
      OW_genes = OW_genes,
      PW_genes = PW_genes,
      mean_surv = mean_surv,
      fecundity_genes = fecundity_genes
    )
  }


#' @title fly method to calculate phenotype
#'
#' @param object An object
#' @param univar_matrix 1D matrix from rcspline.
#' @param bivar_matrix 2D matrix from rcspline.
#' @param logis_k Value for the k parameter for the inverse-logit function.
#' @param logis_x0 Value for the x0 parameter for the inverse-logit function.
#' @returns The same object but modified to have a calculated phenotype.
#' @rdname calculate_phenotype
#' @importFrom rcspline logis spline_2d spline_1d
#' @export
setGeneric("calculate_phenotype", function(object,
                                           univar_matrix,
                                           bivar_matrix,
                                           logis_k,
                                           logis_x0) {
  standardGeneric("calculate_phenotype")
})

#' @title fly method to calculate phenotype
#'
#' @param object An object
#' @param univar_matrix 1D matrix from rcspline.
#' @param bivar_matrix 2D matrix from rcspline.
#' @param logis_k Value for the k parameter for the inverse-logit function.
#' @param logis_x0 Value for the x0 parameter for the inverse-logit function.
#' @returns The same object but modified to have a calculated phenotype.
#' @rdname calculate_phenotype
#' @importFrom rcspline logis spline_2d spline_1d
#' @export
setMethod(
  "calculate_phenotype",
  "fly",
  function(object,
           univar_matrix,
           bivar_matrix,
           logis_k,
           logis_x0) {
    # PE
    object@PE <-
      logis(
        x = spline_2d(object@PE_genes, bivar_matrix),
        max = abs(object@min_PE) + abs(object@max_PE),
        k = logis_k,
        x_0 = logis_x0
      ) + object@min_PE
    # Performing the calculation beforehand instead of inside the
    # ifelse tests.
    # Assigning the values directly to the final object to make it more
    # memory efficient
    object@GW <-
      spline_1d(object@GW_genes, univar_matrix)
    object@OW <-
      spline_1d(object@OW_genes, univar_matrix)
    object@PW <-
      spline_1d(object@PW_genes, univar_matrix)
    # Performing the same truncation as in the original simulation
    object@GW <- ifelse(object@GW < 0,
      0,
      ifelse(object@GW > 1000,
        1000,
        object@GW
      )
    )
    object@OW <- ifelse(object@OW < 0,
      0,
      ifelse(object@OW > 1000,
        1000,
        object@OW
      )
    )
    object@PW <- ifelse(object@PW < 0,
      0,
      ifelse(object@PW > 1000,
        1000,
        object@PW
      )
    )
    weight_sum <- object@GW + object@OW + object@PW
    # Testing for NaN because in some situations some division by zero
    ## can happen. Simply replacing with 0 seems to work
    object@GW <-
      ifelse(is.nan(object@GW / weight_sum), 0, object@GW / weight_sum)
    object@OW <-
      ifelse(is.nan(object@OW / weight_sum), 0, object@OW / weight_sum)
    object@PW <-
      ifelse(is.nan(object@PW / weight_sum), 0, object@PW / weight_sum)
    # egg_laying vector
    # 800 is used as a scale for the maximum number of eggs at the start
    object@egg_laying <-
      800 * skew_normal_pdf(
        object@x_PE,
        object@fecundity_genes[1],
        object@fecundity_genes[2],
        object@fecundity_genes[3]
      )
    return(object)
  }
)

#' @title fly method to retrieve fecundity with day
#'
#' @param object An object
#' @param temperature Temperature for which to retrieve fecundity.
#' @returns Fecundity at the provided temperature.
#' @rdname get_fecundity
#' @export
setGeneric("get_fecundity", function(object,
                                     temperature) {
  standardGeneric("get_fecundity")
})

#' @title fly method to retrieve fecundity with day
#'
#' @param object An object
#' @param temperature Temperature for which to retrieve fecundity.
#' @returns Fecundity at the provided temperature.
#' @rdname get_fecundity
#' @export
setMethod(
  "get_fecundity",
  "fly",
  function(object,
           temperature) {
    output <- numeric()
    for (i in seq_along(temperature)) {
      output <- c(
        output,
        object@egg_laying[
          which.min(abs(object@x_PE - temperature[i]))
        ]
      )
    }
    return(output)
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

#' @title fly method to retrieve prediction with two temperatures
#'
#' @param object An object
#' @param first_temp Value of temperature for first fly experience.
#' @param second_temp Value of temperature for second fly experience.
#' @returns The prediction from the reaction norm.
#' @rdname get_prediction
#' @export
setMethod(
  "get_prediction",
  "fly",
  function(object, first_temp, second_temp) {
    if (length(first_temp) != length(second_temp)) {
      stop("The two temperature objects don't have the same length.")
    } else {
      output <- numeric()
      for (i in seq_along(first_temp)) {
        output <- c(
          output,
          object@PE[
            which.min(abs(object@x_PE - second_temp[i])),
            which.min(abs(object@x_PE - first_temp[i]))
          ]
        )
      }
      return(output)
    }
  }
)

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
#' @returns The survival according to temperature, prediction and phenotype.
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
#' @returns The survival according to temperature, prediction and phenotype.
#' @rdname calculate_survival
#' @importFrom stats rnorm
#' @export
setMethod(
  "calculate_survival",
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
    offspring_difference <-
      abs(offspring_experience - object@mean_surv)
    # in a previous version, I used a different "x_" vector for weights
    # but it has the same length as x_PE with only half the max. This
    # should be equivalent.
    weights_indices <- numeric()
    for (i in length(offspring_difference)) {
      weights_indices <- c(
        weights_indices,
        which.min(abs((object@x_PE / 2) - offspring_difference[i]))
      )
    }
    midpoint <-
      object@mean_surv * overall_GW +
      (1 - overall_GW) *
        (
          object@mean_surv * object@GW[weights_indices] +
            offspring_experience * object@OW[weights_indices] +
            received_PE * object@PW[weights_indices]
        )
    return(c(midpoint, survival(temperature, midpoint)))
  }
)


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
setMethod(
  "get_weight",
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
    offspring_difference <-
      abs(offspring_experience - object@mean_surv)
    # in a previous version, I used a different "x_" vector for weights
    # but it has the same length as x_PE with only half the max. This
    # should be equivalent.
    weights_indices <- numeric()
    for (i in seq_along(offspring_difference)) {
      weights_indices <- c(
        weights_indices,
        which.min(abs((object@x_PE / 2) - offspring_difference[i]))
      )
    }
    if (weight == "GW") {
      return(object@GW[weights_indices])
    } else if (weight == "OW") {
      return(object@OW[weights_indices])
    } else if (weight == "PW") {
      return(object@PW[weights_indices])
    }
  }
)
