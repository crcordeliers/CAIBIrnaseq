# Extract gene expression and sample annotations

This function extracts expression data for a given list of genes from a
SummarizedExperiment object and merges it with sample metadata.

## Usage

``` r
get_exp_df(exp_data, genes, assay = "norm")
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object.

- genes:

  A character vector of gene names or identifiers to extract.

- assay:

  A character string specifying which assay to extract. Defaults to
  \`"norm"\`.

## Value

A data frame with expression values and sample annotations.
