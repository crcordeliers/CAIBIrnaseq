# Score PROGENy Pathways

Calculate PROGENy pathway activity scores based on normalized gene
expression data.

## Usage

``` r
score_progeny(exp_data, species = "Homo sapiens", min_genes = 100)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing normalized gene
  expression data in the \`assays(exp_data)\$norm\` slot.

- species:

  A character string specifying the species, either \`"Homo sapiens"\`
  (default) or \`"Mus musculus"\`.

- min_genes:

  An integer specifying the minimum number of gene symbols in
  \`exp_data\` that must match PROGENy's model genes before scoring
  proceeds. Default is \`100\`.

## Value

A data frame containing PROGENy pathway activity scores for each sample.

## Details

The function uses the normalized gene expression data from the specified
assay of the \`SummarizedExperiment\` object to compute PROGENy scores.
PROGENy scores represent pathway activity and are computed based on
predefined pathway models for the specified species.

PROGENy's models are keyed on gene symbols (HGNC for human, MGI for
mouse), so \`exp_data\` must have gene symbols as rownames of the
\`norm\` assay, not Ensembl (or other) gene IDs. Use \`rebase_gexp()\`
beforehand to convert if needed.
