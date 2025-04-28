#' Transcripts Per Million (TPM) Calculation
#'
#'
#' This function computes the Transcripts Per Million (TPM) normalization for gene expression data.
#' The counts are normalized by the length of each gene and then scaled by the total counts per sample.
#'
#' @param counts A matrix of gene expression counts (rows = genes, columns = samples).
#' @param gene_lengths A named vector of gene lengths. The names must match the row names of `counts`.
#'
#' @returns A matrix of TPM values (rows = genes, columns = samples).
#' @export
#'
transcripts_per_million <- function(counts, gene_lengths) {
  # Check if all gene lengths are provided for the genes in counts
  if(!all(names(gene_lengths) %in% rownames(counts))) {
    stop("All gene lengths need to be provided for the genes in the counts data.")
  }

  # Subset counts to include only genes with provided lengths
  counts <- counts[names(gene_lengths), ]

  # Calculate RPK (Reads Per Kilobase)
  rpk <- counts / gene_lengths

  # Normalize by column (sample) sum and scale to 1 million
  tpm_mat <- t(t(rpk) * 1e6 / colSums(rpk))

  return(tpm_mat)
}
