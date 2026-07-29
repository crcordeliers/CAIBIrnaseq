# Transcripts Per Million (TPM) Calculation

This function computes the Transcripts Per Million (TPM) normalization
for gene expression data. The counts are normalized by the length of
each gene and then scaled by the total counts per sample.

## Usage

``` r
transcripts_per_million(counts, gene_lengths)
```

## Arguments

- counts:

  A matrix of gene expression counts (rows = genes, columns = samples).

- gene_lengths:

  A named vector of gene lengths. The names must match the row names of
  \`counts\`.

## Value

A matrix of TPM values (rows = genes, columns = samples).
