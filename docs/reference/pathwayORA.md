# Pathway Over-Representation Analysis (ORA)

Performs an over-representation analysis (ORA) using Fisher's exact test
on a list of differentially expressed genes against a set of predefined
pathways.

## Usage

``` r
pathwayORA(
  diffexp_result,
  pathways,
  direction = "up",
  id_col = "gene_symbol",
  pcutoff = 0.05
)
```

## Arguments

- diffexp_result:

  A data frame of differential expression results, where row names
  correspond to gene identifiers.

- pathways:

  A data frame with at least two columns: one indicating pathway names
  (e.g., 'pathway') and one with gene identifiers.

- direction:

  either "up" or "down". ORA takes only "up" or "down"-regulated genes
  to do it's overrepresentation analysis.

- id_col:

  Character. The column name in \`pathways\` that matches gene
  identifiers in \`diffexp_result\`.

- pcutoff:

  adjusted p-value to be used to filter differentially expressed genes.

## Value

A data frame of enriched pathways with columns:

- Pathway:

  Name of the enriched pathway

- PValue:

  Raw p-value from Fisher's exact test

- PAdj:

  Adjusted p-value (Benjamini-Hochberg)

- GeneRatio:

  Proportion of input genes found in the pathway

- BgRatio:

  Proportion of background genes found in the pathway

- Genes:

  Comma-separated list of matched genes
