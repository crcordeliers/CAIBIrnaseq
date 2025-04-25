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
estimateTPM <- function(exp, gene_lengths) {

  # --- VALIDATION CHECKS ---

  # Check if exp is a matrix or data.frame
  if (!is.matrix(exp) && !is.data.frame(exp)) {
    stop("Error: `exp` must be a matrix or data frame of numeric counts.")
  }

  # Convert data frame to matrix if needed
  if (is.data.frame(exp)) {
    exp <- as.matrix(exp)
  }

  # Check that all values in exp are numeric
  if (!is.numeric(exp)) {
    stop("Error: `exp` must contain numeric values only.")
  }

  # Check gene_lengths
  if (!is.numeric(gene_lengths)) {
    stop("Error: `gene_lengths` must be a numeric vector.")
  }

  if (length(gene_lengths) != nrow(exp)) {
    stop("Error: Length of `gene_lengths` (", length(gene_lengths),
         ") must match the number of rows in `exp` (", nrow(exp), ").")
  }

  if (any(gene_lengths <= 0)) {
    stop("Error: All values in `gene_lengths` must be positive.")
  }

  # --- TPM CALCULATION ---
  x <- exp / gene_lengths
  tpms <- t(t(x) * 1e6 / colSums(x))

  return(tpms)
}
