# Prepare Gene Expression Data for Heatmap Visualization

This function extracts and formats gene expression data from a
\`SummarizedExperiment\` object for heatmap plotting. It returns a
long-format data frame along with identifiers for samples and features.

## Usage

``` r
prep_exp_hm(expData, genes, assay = "norm", gene_name = "gene_name")
```

## Arguments

- expData:

  A \`SummarizedExperiment\` object.

- genes:

  Character vector of gene names or IDs to extract.

- assay:

  Character. The name of the assay to use (default is "norm").

- gene_name:

  Character. The column name in \`rowData\` containing gene names.
  Default is "gene_name".

## Value

A list containing:

- table:

  A data frame with expression values and sample annotations.

- colv:

  Column variable to use as the x-axis (typically sample IDs).

- rowv:

  Row variable(s) corresponding to the genes.
