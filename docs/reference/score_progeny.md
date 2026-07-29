# Score PROGENy Pathways

Calculate PROGENy pathway activity scores based on normalized gene
expression data.

## Usage

``` r
score_progeny(exp_data, species = "Homo sapiens")
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing normalized gene
  expression data in the \`assays(exp_data)\$norm\` slot.

- species:

  A character string specifying the species, either \`"Homo sapiens"\`
  (default) or \`"Mus musculus"\`.

## Value

A data frame containing PROGENy pathway activity scores for each sample.

## Details

The function uses the normalized gene expression data from the specified
assay of the \`SummarizedExperiment\` object to compute PROGENy scores.
PROGENy scores represent pathway activity and are computed based on
predefined pathway models for the specified species.
