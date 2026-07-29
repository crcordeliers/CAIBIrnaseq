# Identify Highly Variable Genes

This function selects the most highly variable genes from a normalized
gene expression matrix or a \`SummarizedExperiment\` object.

## Usage

``` r
highly_variable_genes(gexp, n_hvg = 2000)
```

## Arguments

- gexp:

  A normalized gene expression matrix (genes in rows, samples in
  columns), or a \`SummarizedExperiment\` object containing normalized
  counts in the \`"norm"\` assay.

- n_hvg:

  An integer specifying the number of highly variable genes to select.
  Default is \`2000\`.

## Value

A character vector containing the names of the top \`n_hvg\` highly
variable genes.

## Details

The function computes the robust coefficient of variation (CV) for each
gene and selects the top \`n_hvg\` genes with the highest CV values. If
the input is a \`SummarizedExperiment\` object, the \`"norm"\` assay is
automatically extracted for processing.
