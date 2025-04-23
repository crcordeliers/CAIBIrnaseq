#' Hierarchical Clustering with Optional PCA Dimensionality Reduction
#'
#' This function performs hierarchical clustering on the input data matrix,
#' optionally reducing dimensionality with PCA first. Supports different distance
#' and linkage methods.
#'
#' @param data A numeric matrix or data frame (features x samples). Rows are features (e.g., genes), columns are samples.
#' @param k Integer. The number of clusters to cut the hierarchical tree into.
#' @param pca Logical. If TRUE (default), perform PCA before clustering.
#' @param n_pcs Integer. Number of principal components to use if `pca = TRUE`. Default is 10.
#' @param dist_method Distance method to use: "euclidean" (default), "pearson", "spearman", etc.
#' @param hc_method Linkage method for hierarchical clustering. Default is "complete".
#'
#' @return A named integer vector with cluster assignments for each sample.
#' @export
#'
#'@importFrom stats as.dist cor cutree dist hclust prcomp
#'
cluster_k_hc <- function(data, k, pca = TRUE, n_pcs = 10,
                         dist_method = "euclidean", hc_method = "complete") {
  if (pca) {
    message("-- Reducing dimensionality using PCA")
    pca_res <- stats::prcomp(t(data), center = TRUE, scale. = TRUE)
    data <- t(pca_res$x[, 1:n_pcs, drop = FALSE])
  }

  if (dist_method %in% c("pearson", "spearman")) {
    dist_mat <- stats::as.dist(1 - stats::cor(data, method = dist_method))
  } else {
    dist_mat <- stats::dist(t(data), method = dist_method)
  }

  message("-- Performing hierarchical clustering into ", k, " groups")
  hc <- stats::hclust(dist_mat, method = hc_method)
  clust_vec <- stats::cutree(hc, k = k)
  return(clust_vec)
}
