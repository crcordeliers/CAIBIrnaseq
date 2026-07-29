# Normalize Gene Expression Data

This function normalizes gene expression counts in a
\`SummarizedExperiment\` object using either \`rlog\` or \`vst\`
normalization, depending on the number of samples.

## Usage

``` r
normalize_gexp(exp_data)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing the gene expression data
  with a \`counts\` assay.

## Value

A \`SummarizedExperiment\` object with an additional \`norm\` assay
containing the normalized expression data.

## Details

The function applies one of two normalization methods: - \`rlog\`
(regularized log transformation) is used for datasets with fewer than 30
samples. - \`vst\` (variance-stabilizing transformation) is used for
datasets with 30 or more samples.

The normalized data is stored in a new assay named \`"norm"\`.
