# Hierarchical Clustering with Optional PCA Dimensionality Reduction

This function performs hierarchical clustering on the input data matrix,
optionally reducing dimensionality with PCA first. Supports different
distance and linkage methods.

## Usage

``` r
cluster_k_hc(
  data,
  k,
  pca = TRUE,
  n_pcs = 10,
  dist_method = "euclidean",
  hc_method = "complete"
)
```

## Arguments

- data:

  A numeric matrix or data frame (features x samples). Rows are features
  (e.g., genes), columns are samples.

- k:

  Integer. The number of clusters to cut the hierarchical tree into.
  Must be a positive integer.

- pca:

  Logical. If TRUE (default), perform PCA before clustering.

- n_pcs:

  Integer. Number of principal components to use if \`pca = TRUE\`.
  Default is 10. Must be a positive integer.

- dist_method:

  Distance method to use: "euclidean" (default), "pearson", "spearman",
  etc. Should be one of the supported methods.

- hc_method:

  Linkage method for hierarchical clustering. Default is "complete".
  Should be one of the supported methods.

## Value

A named integer vector with cluster assignments for each sample.
