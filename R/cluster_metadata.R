#' Cluster Samples Based on Metadata Scores (e.g., Pathways)
#'
#' This function clusters samples using hierarchical clustering on metadata scores
#' (e.g., pathway activity, cell type scores, etc.) stored in the `metadata()`
#' slot of a `SummarizedExperiment` object.
#'
#' @param exp_data A `SummarizedExperiment` object.
#' @param k Integer. Number of clusters to assign.
#' @param metadata_name Character. Name of the metadata slot element to use (e.g., "pathway_scores").
#' @param pca Logical. Whether to reduce dimensionality using PCA before clustering. Default is TRUE.
#' @param n_pcs Integer. Number of principal components to use if PCA is enabled. Default is 10.
#' @param features Optional character vector of row names (e.g., pathways) to use for clustering. If NULL, use all.
#' @param dist_method Distance metric for clustering ("euclidean", "pearson", "spearman", etc.). Default is "euclidean".
#' @param hc_method Linkage method for hierarchical clustering. Default is "complete".
#'
#' @return A `SummarizedExperiment` object with cluster assignments added to `colData()` as `path_cluster`.
#' @export
#'
#' @importFrom SummarizedExperiment colData
#'
cluster_metadata <- function(exp_data, k,
                             metadata_name = "pathway_scores",
                             pca = TRUE,
                             n_pcs = 10,
                             features = NULL,
                             dist_method = "euclidean",
                             hc_method = "complete") {
  scores <- metadata(exp_data)[[metadata_name]]

  if (is.null(scores)) {
    stop(paste0("No scores found in metadata(exp_data)[['", metadata_name, "']]"))
  }

  if (is.null(features)) {
    message("Clustering ", metadata_name, " based on all ", nrow(scores), " scores")
    features <- rownames(scores)
  } else {
    message("Clustering ", metadata_name, " based on provided features")
  }

  clust_res <- cluster_k_hc(scores[features, ], k = k, pca = pca, n_pcs = n_pcs,
                            dist_method = dist_method, hc_method = hc_method)

  sample_annot <- SummarizedExperiment::colData(exp_data)
  sample_annot[["path_cluster"]] <- factor(clust_res)
  SummarizedExperiment::colData(exp_data) <- sample_annot

  return(exp_data)
}
