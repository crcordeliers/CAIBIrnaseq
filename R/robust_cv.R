#' Robust Coefficient of Variation
#'
#' This function computes the robust coefficient of variation for a given vector.
#' The robust coefficient of variation is calculated as the median absolute deviation
#' (MAD) divided by the median of the data, providing a measure of variability
#' that is less sensitive to outliers.
#'
#' @param v A numeric vector of values for which the robust coefficient of variation
#' is calculated.
#'
#' @returns A numeric value representing the robust coefficient of variation for the input vector.
#' @export
#'
#' @importFrom stats mad median
#'
#' @examples
#' robust_cv(c(1, 2, 3, 4, 5))
#' robust_cv(c(1, 1, 1, 1, 1))
#'
robust_cv <- function (v) {
  stats::mad(v, na.rm = TRUE) / stats::median(v, na.rm = TRUE)
}
