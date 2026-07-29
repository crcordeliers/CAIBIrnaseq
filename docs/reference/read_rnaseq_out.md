# Read RNA-seq Output from nf-core/rnaseq

This function loads and formats a \`SummarizedExperiment\` object from
the nf-core/rnaseq output directory (e.g., generated using STAR +
Salmon). It sets up TPM, log-transformed TPM, cleans the sample
annotations, and adds basic QC metadata (number of counts and detected
features).

## Usage

``` r
read_rnaseq_out(DATA_PATH)
```

## Arguments

- DATA_PATH:

  A character string pointing to the path containing nf-core/rnaseq
  results.

## Value

A \`SummarizedExperiment\` object with assays and annotations formatted.
