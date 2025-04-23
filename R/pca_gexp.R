#' Perform PCA on Gene Expression Data
#'
#' This function performs principal component analysis (PCA) on gene expression data in a `SummarizedExperiment` object.
#'
#' @param exp_data A `SummarizedExperiment` object containing the gene expression data.
#' @param assay A character string specifying the assay to use for PCA. Default is `"norm"`.
#' @param filter Logical. If `TRUE`, PCA is performed on the top `n_hvg` highly variable genes. If `FALSE`, PCA is performed on all genes. Default is `TRUE`.
#' @param n_hvg An integer specifying the number of highly variable genes to use if `filter = TRUE`. Default is 2000.
#' @param center Logical. If `TRUE`, the variables are centered before PCA. Default is `TRUE`.
#' @param scale Logical. If `TRUE`, the variables are scaled to unit variance before PCA. Default is `TRUE`.
#'
#' @return A `prcomp` object containing the PCA results.
#'
#' @details
#' If `filter = TRUE`, the function identifies the top `n_hvg` highly variable genes using a robust coefficient of variation and performs PCA on these genes. Otherwise, PCA is performed on all genes in the selected assay.
#'
#' The `prcomp` function is used for PCA, with options to center and scale the data before analysis.
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom stats prcomp
#'
#' @export
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
