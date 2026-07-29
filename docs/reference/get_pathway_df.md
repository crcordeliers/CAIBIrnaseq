# Extract pathway scores and metadata from a SummarizedExperiment

This function extracts pathway activity scores for specified pathways
and merges them with the sample annotations from a SummarizedExperiment
object.

## Usage

``` r
get_pathway_df(exp_data, pathway_scores, pathways)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing sample metadata in
  \`colData\`.

- pathway_scores:

  A matrix or data frame with pathway scores (rows = pathways, columns =
  samples).

- pathways:

  A character vector of pathway names or identifiers to extract.

## Value

A data frame containing the selected pathway scores and sample
annotations.
