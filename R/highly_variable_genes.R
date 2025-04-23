#' Identify Highly Variable Genes
#'
#' This function selects the most highly variable genes from a normalized gene expression matrix or a `SummarizedExperiment` object.
#'
#' @param gexp A normalized gene expression matrix (genes in rows, samples in columns), or a `SummarizedExperiment` object containing normalized counts in the `"norm"` assay.
#' @param n_hvg An integer specifying the number of highly variable genes to select. Default is `2000`.
#'
#' @return A character vector containing the names of the top `n_hvg` highly variable genes.
#'
#' @details
#' The function computes the robust coefficient of variation (CV) for each gene and selects the top `n_hvg` genes with the highest CV values.
#' If the input is a `SummarizedExperiment` object, the `"norm"` assay is automatically extracted for processing.
#'
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
