#' Cluster Samples Based on Gene Expression
#'
#' This function performs hierarchical clustering on samples using normalized gene expression data.
#'
#' @param exp_data A `SummarizedExperiment` object containing the normalized expression data matrix.
#' @param k An integer specifying the number of clusters to generate.
#' @param genes A character vector of gene names to be used for clustering. If `NULL`, the top 2000 highly variable genes are selected automatically.
#' @param pca Logical. If `TRUE`, principal component analysis (PCA) is performed to reduce dimensionality before clustering. Default is `TRUE`.
#' @param n_pcs An integer specifying the number of principal components to retain if `pca = TRUE`. Default is 10.
#' @param dist_method A character string specifying the distance metric to use. Can be one of `"euclidean"`, `"manhattan"`, `"pearson"`, or `"spearman"`. Default is `"euclidean"`.
#' @param hc_method A character string specifying the agglomeration method for hierarchical clustering. Default is `"complete"`.
#'
#' @return A `SummarizedExperiment` object with an updated column in `colData`, named `"exp_cluster"`, containing the cluster assignments as a factor.
#'
#' @details
#' If no specific genes are provided, the function automatically selects the top 2000 highly variable genes for clustering.
#' Optionally, PCA can be applied to reduce dimensionality, which can be useful for datasets with a large number of genes.
#'
#' @importFrom SummarizedExperiment assays colData
#'
#' @export
cluster_exp <- function(exp_data, k, genes = NULL,
                        pca = TRUE,
                        n_pcs = 10,
                        dist_method = "euclidean",
                        hc_method = "complete") {

  # Check if exp_data is a SummarizedExperiment object
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("The 'exp_data' must be a 'SummarizedExperiment' object.")
  }

  # Check if k is a positive integer
  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != round(k)) {
    stop("The 'k' parameter must be a positive integer.")
  }

  # Check if dist_method is valid
  valid_dist_methods <- c("euclidean", "manhattan", "pearson", "spearman")
  if (!(dist_method %in% valid_dist_methods)) {
    stop(paste("The distance method must be one of the following:", paste(valid_dist_methods, collapse = ", ")))
  }

  # Check if hc_method is valid
  valid_hc_methods <- c("complete", "single", "average", "ward.D", "ward.D2", "mcquitty", "median", "centroid")
  if (!(hc_method %in% valid_hc_methods)) {
    stop(paste("The hierarchical clustering method must be one of the following:", paste(valid_hc_methods, collapse = ", ")))
  }

  # If genes is NULL, select top 2000 highly variable genes
  if (is.null(genes)) {
    message("Clustering based on the top 2000 highly variable genes.")
    genes <- highly_variable_genes(exp_data)
  }

  # Check if the genes are valid (present in the data)
  if (!all(genes %in% rownames(exp_data))) {
    stop("Some specified genes are not present in the expression data.")
  }

  # Extract normalized expression data for the selected genes
  gexp <- SummarizedExperiment::assays(exp_data)[["norm"]][genes, ]

  # Perform hierarchical clustering based on the expression data
  clust_res <- cluster_k_hc(gexp,
                            k = k,
                            pca = pca,
                            n_pcs = n_pcs,
                            dist_method = dist_method,
                            hc_method = hc_method)

  # Add clustering results to colData
  sample_annot <- SummarizedExperiment::colData(exp_data)
  sample_annot[["exp_cluster"]] <- factor(clust_res)
  SummarizedExperiment::colData(exp_data) <- sample_annot

  return(exp_data)
}
