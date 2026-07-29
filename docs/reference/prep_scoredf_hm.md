# Prepare Pathway Scores Data Frame for Heatmap

This function prepares a data frame of pathway scores, transposed and
ready for heatmap visualization. It filters and selects a subset of
pathways, either by name or index.

## Usage

``` r
prep_scoredf_hm(pathway_scores, pathways = 1:20)
```

## Arguments

- pathway_scores:

  A matrix or data frame with pathway scores (rows = pathways, columns =
  samples).

- pathways:

  Either a numeric vector of row indices or a character vector of
  pathway names to include. Default is 1:20.

## Value

A list with the following elements:

- table:

  A data frame in wide format (samples in rows, pathways in columns)
  with \`sample_id\` as row identifier.

- colv:

  The column name to use as sample ID (\`"sample_id"\`).

- rowv:

  A character vector of selected pathway names.
