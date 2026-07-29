# Score Pathways

Computes pathway scores for a given expression dataset using specified
scoring methods.

## Usage

``` r
score_pathways(
  exp_data,
  pathways,
  scoring_method = "gsva",
  min_genes = 100,
  verbose = TRUE
)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing normalized expression
  data in the \`assays(exp_data)\$norm\` slot.

- pathways:

  A data frame with pathway definitions, containing at least two
  columns: \`pathway\` (pathway name) and either \`gene_id\` (Ensembl
  IDs) or \`gene_symbol\` (gene symbols).

- scoring_method:

  A character string specifying the scoring method to use. Options are
  \`"gsva"\`, \`"ssgsea"\`, \`"plage"\`, or \`"zscore"\`. Default is
  \`"gsva"\`.

- min_genes:

  An integer specifying the minimum number of genes required for a
  pathway to be considered. Default is \`100\`. Lower only if using
  targetted panel.

- verbose:

  Logical; if \`TRUE\`, prints progress messages during computation.
  Default is \`TRUE\`.

## Value

A data frame of pathway scores, with pathways as row names and samples
as columns. Pathways are sorted by their total score variation.

## Details

The function identifies the gene annotation used in the expression
matrix (\`gene_id\` or \`gene_symbol\`) by matching row names of
\`assays(exp_data)\$norm\` to the \`pathways\` data frame. It then
splits the pathways into gene sets and scores them using the specified
method from the \`GSVA\` package.

The available scoring methods are:

- \`gsva\`:

  Gene Set Variation Analysis.

- \`ssgsea\`:

  Single-sample Gene Set Enrichment Analysis.

- \`plage\`:

  Pathway Level Analysis of Gene Expression.

- \`zscore\`:

  Z-score normalization.

Pathway scores are sorted by the sum of absolute scores across samples,
prioritizing pathways with the highest variation.
