# Prepare Pathway Scores for Heatmap Visualization

This function prepares a tidy table of pathway scores joined with sample
annotations from a \`SummarizedExperiment\` object. It returns a list
suitable for heatmap plotting.

## Usage

``` r
prep_scores_hm(exp_data, pathway_scores, pathways = 1:20)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing sample annotations in
  \`colData\`.

- pathway_scores:

  A matrix or data frame with pathway scores (rows = pathways, columns =
  samples).

- pathways:

  Either a numeric vector indicating indices of pathways to include, or
  a character vector of pathway names. Default is 1:20.

## Value

A list containing:

- table:

  A data frame in wide format with pathway scores and sample metadata.

- colv:

  Column to use as identifier for samples (usually \`"sample_id"\`).

- rowv:

  A character vector of selected pathways.
