# fly.R

#' Fly class
#' @description
#' Class to keep individual values together inside a single "fly".
#' Will be created with the specific values for each genomic trait and the
#' phenotype will be calculated with internal functions
#'
#' @slot replicate The simulation replicate the fly was in.
#' @slot population Population in which the fly was in.
#' @slot PE_genes Vector of 25 values that define the PE reaction norm.
#' @slot GW_genes Vector of 5 vales that define the GW reaction norm.
#' @slot OW_genes Vector of 5 vales that define the OW reaction norm.
#' @slot PW_genes Vector of 5 vales that define the PW reaction norm.
#' @slot mean_surv Genetic value for the survival midpoint.
#' @slot fecundity_genes Vector of 3 values that define the mean, sd and skew for fecundity.
#' @slot PE PE reaction norm.
#' @slot GW GW reaction norm.
#' @slot OW OW reaction norm.
#' @slot PW PW reaction norm.
#' @slot egg_laying Reaction for egg laying.
setClass(
  "fly",
  slots = list(
    replicate = "integer",
    population = "integer",
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
    replicate = 0L,
    population = 0L,
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
    fecundity_genes = c(25, 15,-5),
    # Below are phenotypes
    PE = matrix(0, 100, 100),
    GW = rep(0, 100),
    OW = rep(0, 100),
    PW = rep(0, 100),
    egg_laying = rep(0, 100)
  )
)

# Constructor method for fly
fly <-
  function(replicate,
           population,
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
#' @param max_PE Numeric value that represents the maximum PE.
#' @param logis_k Value for the k parameter for the inverse-logit function.
#' @param logis_x0 Value for the x0 parameter for the inverse-logit function.
#' @param precision Value defining how precise the calculations should be.
#' @rdname calculate_phenotype
setGeneric("calculate_phenotype",  function(object,
                                            univar_matrix,
                                            bivar_matrix,
                                            max_PE,
                                            logis_k,
                                            logis_x0,
                                            precision)
  standardGeneric("calculate_phenotype"))

setMethod("calculate_phenotype",
          "fly",
          function(object,
                   univar_matrix,
                   bivar_matrix,
                   max_PE,
                   logis_k,
                   logis_x0,
                   precision) {
            # PE
            object@PE <-
              rcspline::logis(
                x = rcspline::spline_2d(object@PE_genes, bivar_matrix),
                max = max_PE,
                k = logis_k,
                x_0 = logis_x0
              )
            
            # Performing the calculation beforehand instead of inside the ifelse tests.
            # Assigning the values directly to the final object to make it more
            # memory efficient
            object@GW <-
              rcspline::spline_1d(object@GW_genes, univar_matrix)
            object@OW <-
              rcspline::spline_1d(object@OW_genes, univar_matrix)
            object@PW <-
              rcspline::spline_1d(object@PW_genes, univar_matrix)
            
            # Performing the same truncation as in the original simulation
            object@GW <- ifelse(object@GW < 0,
                                0,
                                ifelse(object@GW > 1000,
                                       1000,
                                       object@GW))
            object@OW <- ifelse(object@OW < 0,
                                0,
                                ifelse(object@OW > 1000,
                                       1000,
                                       object@OW))
            object@PW <- ifelse(object@PW < 0,
                                0,
                                ifelse(object@PW > 1000,
                                       1000,
                                       object@PW))
            
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
            object@egg_laying <-
              fecundity(
                seq(0, max_PE, length = precision),
                object@fecundity_genes[1],
                object@fecundity_genes[2],
                object@fecundity_genes[3]
              )
          })
