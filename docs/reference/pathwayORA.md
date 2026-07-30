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

  Either a data frame with at least two columns, one indicating pathway
  names (\`pathway\`) and one with gene identifiers (named by
  \`id_col\`), or a named list of character vectors (one gene-ID vector
  per pathway/signature), e.g. \`list(MySignature = c("GENE1",
  "GENE2"))\`. The list form lets custom gene signatures be tested the
  same way as an MSigDB collection, and ignores \`id_col\`.

- direction:

  either "up" or "down". ORA takes only "up" or "down"-regulated genes
  to do it's overrepresentation analysis.

- id_col:

  Character. The column name in \`pathways\` that matches gene
  identifiers in \`diffexp_result\`. Only used when \`pathways\` is a
  data frame.

- pcutoff:

  adjusted p-value to be used to filter differentially expressed genes.

## Value

A data frame of enriched pathways with columns:

- pathway:

  Name of the enriched pathway

- pval:

  Raw p-value from Fisher's exact test

- padj:

  Adjusted p-value (Benjamini-Hochberg)

- geneRatio:

  Proportion of input genes found in the pathway

- bgRatio:

  Proportion of background genes found in the pathway

- genes:

  Comma-separated list of matched genes

Column names match the naming used by
[`pathwayFGSEA`](https://crcordeliers.github.io/CAIBIrnaseq/reference/pathwayFGSEA.md)
(`pathway`, `pval`, `padj`) so that downstream consumers such as
[`plot_pathway_dotplot`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_pathway_dotplot.md)
can handle both result types.
