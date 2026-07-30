# Pathway Enrichment Using FGSEA

Performs fast gene set enrichment analysis (FGSEA) using the
\`fgseaMultilevel\` method, based on a ranked list of differentially
expressed genes and a collection of pathways.

## Usage

``` r
pathwayFGSEA(diffexp, pathwayCollection, seed = 0)
```

## Arguments

- diffexp:

  A data frame of differential expression results, with row names as
  gene IDs and a \`log2FoldChange\` column used for ranking genes.

- pathwayCollection:

  Either a data frame with at least two columns, one for pathway names
  (\`pathway\`) and one for gene symbols (\`gene_symbol\`), or a named
  list of character vectors (one gene-symbol vector per
  pathway/signature), e.g. \`list(MySignature = c("GENE1", "GENE2"))\`.
  The list form lets custom gene signatures be tested the same way as an
  MSigDB collection.

## Value

A data frame returned by \`fgseaMultilevel()\`, including columns such
as:

- pathway:

  Name of the pathway

- pval:

  P-value of enrichment

- padj:

  Adjusted p-value (FDR)

- ES:

  Enrichment Score

- NES:

  Normalized Enrichment Score

- leadingEdge:

  Vector of leading-edge genes
