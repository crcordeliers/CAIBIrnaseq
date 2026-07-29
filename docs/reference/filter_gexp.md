# Filter Gene Expression Data

This function filters genes from a \`SummarizedExperiment\` object based
on minimum expression thresholds.

## Usage

``` r
filter_gexp(exp_data, min_nsamp = 1, min_counts = 10)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing the gene expression data
  with a \`counts\` assay.

- min_nsamp:

  An integer specifying the minimum number of samples in which a gene
  must have expression above \`min_counts\` to be retained. Default is
  1.

- min_counts:

  An integer specifying the minimum count threshold a gene must have in
  \`min_nsamp\` samples to be retained. Default is 10.

## Value

A filtered \`SummarizedExperiment\` object containing only the genes
that meet the specified criteria. Two new columns are added to
\`colData\`: - \`ncounts\`: Total counts per sample. - \`nfeats\`:
Number of features (genes) detected per sample.

## Details

This function removes genes that do not meet the specified thresholds
for expression. It adds sample-level metrics (\`ncounts\` and
\`nfeats\`) to the \`colData\` of the \`SummarizedExperiment\` object
for downstream analysis.

The filtering criteria are: - A gene must have expression greater than
or equal to \`min_counts\` in at least \`min_nsamp\` samples.
