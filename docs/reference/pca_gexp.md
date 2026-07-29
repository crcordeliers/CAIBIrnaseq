# Perform PCA on Gene Expression Data

This function performs principal component analysis (PCA) on gene
expression data in a \`SummarizedExperiment\` object.

## Usage

``` r
pca_gexp(
  exp_data,
  assay = "norm",
  filter = TRUE,
  n_hvg = 2000,
  center = TRUE,
  scale = TRUE
)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing the gene expression data.

- assay:

  A character string specifying the assay to use for PCA. Default is
  \`"norm"\`.

- filter:

  Logical. If \`TRUE\`, PCA is performed on the top \`n_hvg\` highly
  variable genes. If \`FALSE\`, PCA is performed on all genes. Default
  is \`TRUE\`.

- n_hvg:

  An integer specifying the number of highly variable genes to use if
  \`filter = TRUE\`. Default is 2000.

- center:

  Logical. If \`TRUE\`, the variables are centered before PCA. Default
  is \`TRUE\`.

- scale:

  Logical. If \`TRUE\`, the variables are scaled to unit variance before
  PCA. Default is \`TRUE\`.

## Value

A \`prcomp\` object containing the PCA results.

## Details

If \`filter = TRUE\`, the function identifies the top \`n_hvg\` highly
variable genes using a robust coefficient of variation and performs PCA
on these genes. Otherwise, PCA is performed on all genes in the selected
assay.

The \`prcomp\` function is used for PCA, with options to center and
scale the data before analysis.
