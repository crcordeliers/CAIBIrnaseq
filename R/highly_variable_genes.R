#' Identify Highly Variable Genes (HVGs)
#'
#' This function identifies the top n highly variable genes (HVGs) from a gene expression matrix.
#' It calculates the robust coefficient of variation for each gene and selects the top `n_hvg` genes with the highest variability.
#'
#' @param gexp A gene expression matrix or a SummarizedExperiment object containing normalized gene expression data.
#' @param n_hvg The number of highly variable genes to return (default is 2000).
#'
#' @returns A character vector of gene names representing the highly variable genes.
#' @export
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom stats na.omit
#' @importFrom magrittr %>%
#'
highly_variable_genes <- function(gexp, n_hvg = 2000) {
  # Check if input is a SummarizedExperiment and extract normalized expression if so
  if("SummarizedExperiment" %in% class(gexp)) {
    if (!"norm" %in% names(SummarizedExperiment::assays(gexp))) {
      stop("Normalized expression data ('norm') not found in the SummarizedExperiment object.")
    }
    gexp <- SummarizedExperiment::assays(gexp)[["norm"]]  # Use the normalized gene expression data
  }

  # Calculate the robust coefficient of variation (CV) for each gene
  gcvs <- apply(gexp, 1, robust_cv)

  # Get the top n highly variable genes based on the CV
  gkeep <- gcvs %>%
    sort(decreasing = TRUE) %>%
    .[1:n_hvg] %>%
    names() %>%
    na.omit()  # Remove NA values

  return(gkeep)
}
