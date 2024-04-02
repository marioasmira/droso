# constant_area.R


#' Class to hold the initial area as a constant
#'
#' @slot lower_bound Value for the lower bound of integration.
#' @slot upper_bound Value for the upper bound of integration.
#' @slot initial_area Value of area that will be kept constant
#'
#' @export
setClass(
  "constant_area",
  slots = list(
    lower_bound = "numeric",
    upper_bound = "numeric",
    initial_area = "numeric"
  ),
  prototype = list(
    lower_bound = -10,
    upper_bound = 40,
    initial_area = 0
  )
)

#' Constructor for constant_area
#'
#' @param lower_bound Value for the lower bound of integration.
#' @param upper_bound Value for the upper bound of integration.
#' @param func Function to use for calculation
#' @param ... Arguments for the function
#'
#' @importFrom stats integrate
#' @return a constant_area object
#' @export
constant_area <- function(lower_bound,
                          upper_bound,
                          func,
                          ...) {
  f <- function(x) {
    ifelse(
      func(x, ...) < 0,
      0,
      func(x, ...)
    )
  }
  area <- integrate(f, lower_bound, upper_bound)
  new(
    "constant_area",
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    initial_area = unlist(area[1])
  )
}

#' Method to compare areas to initial area
#'
#' @param object The constant_area object to use.
#' @param func Function to use for calculation
#' @param ... Arguments for the function
#'
#' @importFrom stats integrate
#' @return Numeric scale of the comparison betwee initial and provided areas.
#' @rdname get_scale
#' @export
setGeneric(
  "get_scale",
  function(func,
           ...) {
    standardGeneric("get_scale")
  }
)

#' @rdname get_scale
setMethod(
  "get_scale",
  "constant_area",
  definition = function(object,
                        func,
                        ...) {
    f <- function(x) {
      ifelse(
        func(x, ...) < 0,
        0,
        func(x, ...)
      )
    }
    area <- tryCatch(
      integrate(f, object@lower_bound, object@upper_bound),
      error = function(e) {
        print(e)
        for (arg in list(...)) {
          print(paste0("value: ", arg, "."))
        }
      }
    )

    return(object@initial_area / unlist(area[1]))
  }
)
