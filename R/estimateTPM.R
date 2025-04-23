#' Estimate Transcripts Per Million (TPM)
#'
#' This function calculates TPM values from raw expression counts and gene lengths.
#' TPM is a normalization method that accounts for both sequencing depth and gene length,
#' allowing comparison of gene expression levels within and between samples.
#'
#' @param exp A numeric matrix or data frame of raw expression counts. Rows represent genes, columns represent samples.
#' @param gene_lengths A numeric vector of gene lengths (in kilobases or base pairs) corresponding to the rows of \code{exp}.
#'                     Must be the same length as the number of rows in \code{exp}.
#'
#' @return A numeric matrix of TPM-normalized expression values with the same dimensions as \code{exp}.
#'
#' @export
#'
estimateTPM <- function(exp, gene_lengths) {
  x <- exp / gene_lengths
  tpms <- t(t(x) * 1e6 / colSums(x))
  return(tpms)
}
