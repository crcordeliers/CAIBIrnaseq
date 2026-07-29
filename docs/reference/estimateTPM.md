# Estimate Transcripts Per Million (TPM)

This function calculates TPM values from raw expression counts and gene
lengths. TPM is a normalization method that accounts for both sequencing
depth and gene length, allowing comparison of gene expression levels
within and between samples.

## Usage

``` r
estimateTPM(exp, gene_lengths)
```

## Arguments

- exp:

  A numeric matrix or data frame of raw expression counts. Rows
  represent genes, columns represent samples.

- gene_lengths:

  A numeric vector of gene lengths (in kilobases or base pairs)
  corresponding to the rows of `exp`. Must be the same length as the
  number of rows in `exp`.

## Value

A numeric matrix of TPM-normalized expression values with the same
dimensions as `exp`.
