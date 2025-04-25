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

  # ---- VALIDATION CHECKS ----

  # Check that exp_data is a SummarizedExperiment
  if (!is(exp_data, "SummarizedExperiment")) {
    stop("Error: `exp_data` must be a SummarizedExperiment object.")
  }

  # Check that metadata_name exists in metadata
  scores <- S4Vectors::metadata(exp_data)[[metadata_name]]
  if (is.null(scores)) {
    stop(paste0("Error: No scores found in metadata(exp_data)[['", metadata_name, "']]"))
  }

  # Check that k is a positive integer
  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != floor(k)) {
    stop("Error: `k` must be a positive integer.")
  }

  # Check that scores is a matrix or data.frame
  if (!is.matrix(scores) && !is.data.frame(scores)) {
    stop("Error: The metadata `", metadata_name, "` must be a matrix or data frame.")
  }

  # If features are specified, check that they exist in scores
  if (!is.null(features)) {
    missing_features <- setdiff(features, rownames(scores))
    if (length(missing_features) > 0) {
      stop("Error: The following features are not found in metadata rows: ",
           paste(missing_features, collapse = ", "))
    }
    message("Clustering ", metadata_name, " based on provided features (", length(features), ")")
  } else {
    features <- rownames(scores)
    message("Clustering ", metadata_name, " based on all ", length(features), " features")
  }

  # Check n_pcs if PCA is TRUE
  if (pca) {
    if (!is.numeric(n_pcs) || length(n_pcs) != 1 || n_pcs <= 0 || n_pcs != floor(n_pcs)) {
      stop("Error: `n_pcs` must be a positive integer.")
    }
    if (n_pcs > length(features)) {
      warning("Warning: `n_pcs` is greater than the number of selected features. Using all available features.")
      n_pcs <- length(features)
    }
  }

  # Validate dist_method
  valid_dist_methods <- c("euclidean", "pearson", "spearman")
  if (!(dist_method %in% valid_dist_methods)) {
    stop("Error: `dist_method` must be one of: ", paste(valid_dist_methods, collapse = ", "))
  }

  # Validate hc_method
  valid_hc_methods <- c("complete", "single", "average", "ward.D", "ward.D2")
  if (!(hc_method %in% valid_hc_methods)) {
    stop("Error: `hc_method` must be one of: ", paste(valid_hc_methods, collapse = ", "))
  }

  # ---- CLUSTERING ----

  clust_res <- cluster_k_hc(scores[features, ], k = k, pca = pca, n_pcs = n_pcs,
                            dist_method = dist_method, hc_method = hc_method)

  sample_annot <- SummarizedExperiment::colData(exp_data)
  sample_annot[["path_cluster"]] <- factor(clust_res)
  SummarizedExperiment::colData(exp_data) <- sample_annot

  return(exp_data)
}
