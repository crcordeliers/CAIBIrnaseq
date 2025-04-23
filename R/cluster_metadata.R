#' Cluster Metadata in Expression Data
#'
#' This function performs hierarchical clustering on specified metadata scores
#' within the provided expression dataset. It updates the dataset's sample annotations
#' with the clustering results.
#'
#' @param exp_data A SummarizedExperiment object containing expression data and metadata.
#' @param k An integer specifying the number of clusters to generate.
#' @param metadata_name A character string specifying the name of the metadata variable
#'   to use for clustering. Defaults to "pathway_scores".
#' @param pca Logical. If TRUE, principal component analysis (PCA) is applied before clustering. Defaults to TRUE.
#' @param n_pcs An integer specifying the number of principal components to retain if PCA is applied. Defaults to 10.
#' @param features A character vector of feature names to use for clustering. If NULL, all features are used. Defaults to NULL.
#' @param dist_method A character string specifying the distance metric to use for hierarchical clustering. Defaults to "euclidean".
#' @param hc_method A character string specifying the hierarchical clustering linkage method. Defaults to "complete".
#'
#' @return A SummarizedExperiment object with updated sample annotations including the clustering results in the "path_cluster" column.
#'
#' @details If no features are provided, all rows in the specified metadata are used for clustering.
#'   PCA can be applied to reduce dimensionality before clustering.
#'
#' @importFrom SummarizedExperiment colData
#'
#' @export
cluster_metadata <- function(exp_data, k,
                             metadata_name = "pathway_scores",
                             pca = TRUE,
                             n_pcs = 10,
                             features = NULL,
                             dist_method = "euclidean",
                             hc_method = "complete") {
  scores <- S4Vectors::metadata(exp_data)[[metadata_name]]

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
