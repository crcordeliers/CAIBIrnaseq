# Cluster Samples Based on Gene Expression

This function performs hierarchical clustering on samples using
normalized gene expression data.

## Usage

``` r
cluster_exp(
  exp_data,
  k,
  genes = NULL,
  pca = TRUE,
  n_pcs = 10,
  dist_method = "euclidean",
  hc_method = "complete"
)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing the normalized expression
  data matrix.

- k:

  An integer specifying the number of clusters to generate.

- genes:

  A character vector of gene names to be used for clustering. If
  \`NULL\`, the top 2000 highly variable genes are selected
  automatically.

- pca:

  Logical. If \`TRUE\`, principal component analysis (PCA) is performed
  to reduce dimensionality before clustering. Default is \`TRUE\`.

- n_pcs:

  An integer specifying the number of principal components to retain if
  \`pca = TRUE\`. Default is 10.

- dist_method:

  A character string specifying the distance metric to use. Can be one
  of \`"euclidean"\`, \`"manhattan"\`, \`"pearson"\`, or \`"spearman"\`.
  Default is \`"euclidean"\`.

- hc_method:

  A character string specifying the agglomeration method for
  hierarchical clustering. Default is \`"complete"\`.

## Value

A \`SummarizedExperiment\` object with an updated column in \`colData\`,
named \`"exp_cluster"\`, containing the cluster assignments as a factor.

## Details

If no specific genes are provided, the function automatically selects
the top 2000 highly variable genes for clustering. Optionally, PCA can
be applied to reduce dimensionality, which can be useful for datasets
with a large number of genes.
