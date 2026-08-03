# Pathway Analysis

Performs pathway enrichment analysis using either Over-representation
Analysis (ORA) or Gene Set Enrichment Analysis (FGSEA) based on the
provided gene set and pathways.

## Usage

``` r
pathwayAnalysis(
  diffexp,
  pathways,
  method = "ORA",
  id_col = "gene_symbol",
  pcutoff = 0.05,
  ...
)
```

## Arguments

- pathways:

  Either a data frame containing pathway information, with columns for
  pathways and associated genes, or a named list of character vectors
  (one gene-ID vector per pathway/signature), e.g. \`list(MySignature =
  c("GENE1", "GENE2"))\`. The list form lets custom gene signatures be
  tested/scored the same way as an MSigDB collection.

- method:

  A string specifying the pathway analysis method. Options are \`"ORA"\`
  for Over-representation Analysis or \`"FGSEA"\` for Gene Set
  Enrichment Analysis. Default is \`"ORA"\`.

- id_col:

  A string specifying the column in the \`pathways\` data frame that
  contains gene identifiers. Default is \`"gene_symbol"\`. Only used
  when \`pathways\` is a data frame.

- pcutoff:

  A numeric value specifying the adjusted p-value cutoff for significant
  pathways. Default is \`0.05\`.

- gene_set:

  A data frame containing the differential expression results. Rows
  represent genes.

- verbose:

  A logical value indicating whether to display progress messages.
  Default is \`TRUE\`.

## Value

A data frame containing pathway enrichment results, including pathway
names, p-values, adjusted p-values, and other relevant statistics.

## Details

This function allows users to perform pathway enrichment analysis using
either ORA or FGSEA. The \`pathways\` input should either be a data
frame including a column specifying gene identifiers and a column
specifying pathway names, or a named list of gene-ID vectors.

\- \*\*ORA\*\*: Identifies pathways over-represented in the provided
gene set. - \*\*FGSEA\*\*: Identifies pathways using a ranked list of
genes based on differential expression statistics.

The function dynamically calls internal implementations of ORA or FGSEA
based on the \`method\` argument.
