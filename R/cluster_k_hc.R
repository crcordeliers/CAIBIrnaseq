#' Hierarchical Clustering with Optional PCA Dimensionality Reduction
#'
#' This function performs hierarchical clustering on the input data matrix,
#' optionally reducing dimensionality with PCA first. Supports different distance
#' and linkage methods.
#'
#' @param data A numeric matrix or data frame (features x samples). Rows are features (e.g., genes), columns are samples.
#' @param k Integer. The number of clusters to cut the hierarchical tree into. Must be a positive integer.
#' @param pca Logical. If TRUE (default), perform PCA before clustering.
#' @param n_pcs Integer. Number of principal components to use if `pca = TRUE`. Default is 10. Must be a positive integer.
#' @param dist_method Distance method to use: "euclidean" (default), "pearson", "spearman", etc. Should be one of the supported methods.
#' @param hc_method Linkage method for hierarchical clustering. Default is "complete". Should be one of the supported methods.
#'
#' @return A named integer vector with cluster assignments for each sample.
#' @export
#'
#' @importFrom stats as.dist cor cutree dist hclust prcomp
#'
cluster_k_hc <- function(data, k, pca = TRUE, n_pcs = 10,
                         dist_method = "euclidean", hc_method = "complete") {

  # Check if `data` is a numeric matrix or data frame
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Error: `data` must be a numeric matrix or data frame.")
  }

  # Check if the number of clusters `k` is a positive integer
  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != floor(k)) {
    stop("Error: `k` must be a positive integer.")
  }

  # Check if `n_pcs` is a positive integer when PCA is performed
  if (pca) {
    if (!is.numeric(n_pcs) || length(n_pcs) != 1 || n_pcs <= 0 || n_pcs != floor(n_pcs)) {
      stop("Error: `n_pcs` must be a positive integer.")
    }
    if (n_pcs > nrow(data)) {
      warning("Warning: `n_pcs` is greater than the number of features. Using all features.")
      n_pcs <- nrow(data)
    }
  }

  # Check if `dist_method` is valid
  valid_dist_methods <- c("euclidean", "pearson", "spearman")
  if (!(dist_method %in% valid_dist_methods)) {
    stop("Error: `dist_method` must be one of the following: ", paste(valid_dist_methods, collapse = ", "), ".")
  }

  # Check if `hc_method` is valid
  valid_hc_methods <- c("complete", "single", "average", "ward.D", "ward.D2")
  if (!(hc_method %in% valid_hc_methods)) {
    stop("Error: `hc_method` must be one of the following: ", paste(valid_hc_methods, collapse = ", "), ".")
  }

  # Perform PCA if requested
  if (pca) {
    message("-- Reducing dimensionality using PCA")
    pca_res <- stats::prcomp(t(data), center = TRUE, scale. = TRUE)
    data <- t(pca_res$x[, 1:n_pcs, drop = FALSE])
  }

  # Calculate the distance matrix based on the chosen method
  if (dist_method %in% c("pearson", "spearman")) {
    dist_mat <- stats::as.dist(1 - stats::cor(data, method = dist_method))
  } else {
    dist_mat <- stats::dist(t(data), method = dist_method)
  }

  # Perform hierarchical clustering
  message("-- Performing hierarchical clustering into ", k, " groups")
  hc <- stats::hclust(dist_mat, method = hc_method)

  # Cut the tree into `k` clusters
  clust_vec <- stats::cutree(hc, k = k)

  return(clust_vec)
}
