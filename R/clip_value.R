# clip_value.R

#' Clip value to between bounds
#'
#' @param x Value to clip.
#' @param l Lower bound for clipping.
#' @param u Upper bound for clipping.
#'
#' @return Clipped value between l and u.
#' @export
clip_value <- function(
    x,
    l,
    u
    ) {
  ifelse(
    x < l,
    l,
    ifelse(
      x > u,
      u,
      x
    )
  )
}
