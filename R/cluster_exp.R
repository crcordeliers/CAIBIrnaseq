#' Cluster Expression Data Using Hierarchical Clustering
#'
#' This function performs hierarchical clustering on normalized gene expression data
#' from a `SummarizedExperiment` object. Clustering is based on a selected set of genes,
#' optionally using PCA for dimensionality reduction. The resulting cluster assignments
#' are added to the sample metadata (`colData`) as a new column `exp_cluster`.
#'
#' @param exp_data A `SummarizedExperiment` object with normalized expression data stored in the "norm" assay.
#' @param k Integer. The number of clusters to assign.
#' @param genes Optional character vector of gene IDs or symbols to use for clustering. If NULL (default), the top 2000 highly variable genes are selected using `highly_variable_genes()`.
#' @param pca Logical. Whether to reduce dimensionality using PCA before clustering. Default is TRUE.
#' @param n_pcs Integer. Number of principal components to retain if `pca = TRUE`. Default is 10.
#' @param dist_method Character. Distance metric to use ("euclidean", "pearson", "spearman", etc.). Default is "euclidean".
#' @param hc_method Character. Hierarchical clustering linkage method (e.g., "complete", "average", "ward.D"). Default is "complete".
#'
#' @return A `SummarizedExperiment` object with an added `exp_cluster` column in `colData`.
#' @export
#'
#' @importFrom SummarizedExperiment colData
#'
cluster_exp <- function(exp_data, k, genes = NULL,
                        pca = TRUE,
                        n_pcs = 10,
                        dist_method = "euclidean",
                        hc_method = "complete") {
  if (is.null(genes)) {
    message("Clustering based on top 2000 highly variable genes")
    genes <- highly_variable_genes(exp_data)
  }

  gexp <- assays(exp_data)[["norm"]][genes, ]

  clust_res <- cluster_k_hc(gexp,
                            k = k,
                            pca = pca,
                            n_pcs = n_pcs,
                            dist_method = dist_method,
                            hc_method = hc_method)

  sample_annot <- SummarizedExperiment::colData(exp_data)
  sample_annot[["exp_cluster"]] <- factor(clust_res)
  SummarizedExperiment::colData(exp_data) <- sample_annot

  return(exp_data)
}
