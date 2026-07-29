# Cluster Metadata in Expression Data

This function performs hierarchical clustering on specified metadata
scores within the provided expression dataset. It updates the dataset's
sample annotations with the clustering results.

## Usage

``` r
cluster_metadata(
  exp_data,
  k,
  metadata_name = "pathway_scores",
  pca = TRUE,
  n_pcs = 10,
  features = NULL,
  dist_method = "euclidean",
  hc_method = "complete"
)
```

## Arguments

- exp_data:

  A SummarizedExperiment object containing expression data and metadata.

- k:

  An integer specifying the number of clusters to generate.

- metadata_name:

  A character string specifying the name of the metadata variable to use
  for clustering. Defaults to "pathway_scores".

- pca:

  Logical. If TRUE, principal component analysis (PCA) is applied before
  clustering. Defaults to TRUE.

- n_pcs:

  An integer specifying the number of principal components to retain if
  PCA is applied. Defaults to 10.

- features:

  A character vector of feature names to use for clustering. If NULL,
  all features are used. Defaults to NULL.

- dist_method:

  A character string specifying the distance metric to use for
  hierarchical clustering. Defaults to "euclidean".

- hc_method:

  A character string specifying the hierarchical clustering linkage
  method. Defaults to "complete".

## Value

A SummarizedExperiment object with updated sample annotations including
the clustering results in the "path_cluster" column.

## Details

If no features are provided, all rows in the specified metadata are used
for clustering. PCA can be applied to reduce dimensionality before
clustering.
