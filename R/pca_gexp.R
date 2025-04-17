#' Principal Component Analysis (PCA) on Gene Expression Data
#'
#' This function performs PCA (Principal Component Analysis) on the gene expression data.
#' It can optionally filter for highly variable genes (HVGs) before performing PCA.
#'
#' @param exp_data A `SummarizedExperiment` object containing the gene expression data.
#' @param assay A character string specifying the assay to use (default is `"norm"` for normalized data).
#' @param filter Logical value, whether to filter for highly variable genes (default is `TRUE`).
#' @param n_hvg The number of highly variable genes to select if `filter` is `TRUE` (default is 2000).
#' @param center Logical value, whether to center the data before performing PCA (default is `TRUE`).
#' @param scale Logical value, whether to scale the data before performing PCA (default is `TRUE`).
#'
#' @returns A list containing the results of the PCA, including rotation, standard deviation, and scores.
#' @export
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom stats prcomp
#'
pca_gexp <- function(exp_data, assay = "norm", filter = TRUE, n_hvg = 2000, center = TRUE,
                     scale = TRUE) {
  # Extract the gene expression data from the specified assay
  gexp <- SummarizedExperiment::assays(exp_data)[[assay]]

  # Filter for highly variable genes if specified
  if (filter) {
    gkeep <- highly_variable_genes(gexp, n_hvg)
  } else {
    gkeep <- rownames(gexp)
  }

  # Perform PCA on the filtered gene expression data (centered and scaled)
  pca_res <- stats::prcomp(t(gexp[gkeep,]), center = center, scale = scale)

  return(pca_res)
}
