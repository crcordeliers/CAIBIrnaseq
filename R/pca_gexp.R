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
  # Basic checks
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("`exp_data` must be a 'SummarizedExperiment' object.")
  }

  if (!is.character(assay) || length(assay) != 1) {
    stop("`assay` must be a character string of length 1.")
  }

  if (!(assay %in% names(SummarizedExperiment::assays(exp_data)))) {
    stop(paste0("Assay '", assay, "' not found in `exp_data`."))
  }

  if (!is.logical(filter) || length(filter) != 1) {
    stop("`filter` must be a logical value.")
  }

  if (!is.numeric(n_hvg) || length(n_hvg) != 1 || n_hvg <= 0 || n_hvg != as.integer(n_hvg)) {
    stop("`n_hvg` must be a strictly positive integer.")
  }

  if (!is.logical(center) || length(center) != 1) {
    stop("`center` must be a logical value.")
  }

  if (!is.logical(scale) || length(scale) != 1) {
    stop("`scale` must be a logical value.")
  }

  # Extract gene expression matrix
  gexp <- SummarizedExperiment::assays(exp_data)[[assay]]

  if (!is.matrix(gexp)) {
    stop("The extracted assay data must be a matrix.")
  }

  # Gene filtering
  if (filter) {
    gkeep <- highly_variable_genes(gexp, n_hvg)
  } else {
    gkeep <- rownames(gexp)
  }

  if (length(gkeep) == 0) {
    stop("No genes were retained for PCA analysis.")
  }

  # PCA
  pca_res <- stats::prcomp(t(gexp[gkeep,]), center = center, scale. = scale)

  return(pca_res)
}
