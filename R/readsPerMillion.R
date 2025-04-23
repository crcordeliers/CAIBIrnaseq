#' Normalize Counts to Reads Per Million (RPM)
#'
#' This function normalizes raw read counts to Reads Per Million (RPM) across each column (typically samples).
#' It is commonly used for comparing sequencing depth-normalized gene expression across samples.
#'
#' @param data A numeric matrix or data frame of raw counts. Rows are typically genes, and columns are samples.
#' @param factor A numeric value for scaling. Default is \code{10^6} to compute reads per million.
#'
#' @return A matrix of the same dimensions as \code{data}, with RPM-normalized values.
#'
#' @export
#'
readsPerMillion <- function(data, factor = 10^6) {
  apply(data, 2, function(col) col / sum(col) * factor)
}
