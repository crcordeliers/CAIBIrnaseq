# Plot Volcano Plot of Differential Expression Results

Generates a volcano plot based on the results of differential expression
analysis, highlighting upregulated and downregulated genes, with labels
for top significant genes.

## Usage

``` r
plot_exp_volcano(
  diffexp,
  nb = 10,
  title = "Volcano Plot of Differential Expression",
  color_up = "#0072B2",
  color_down = "#D55E00",
  color_ns = "gray80",
  fname = NULL
)
```

## Arguments

- diffexp:

  A data frame containing differential expression results. Must
  include: - \`log2FoldChange\`: The log2 fold change values for each
  gene. - \`padj\`: The adjusted p-value for each gene.

- nb:

  The number of genes that have an annotation

- color_up:

  Color used for upregulated genes (default is "#0072B2").

- color_down:

  Color used for downregulated genes (default is "#D55E00").

- color_ns:

  Color used for non-significant genes (default is "gray80").

## Value

A \`ggplot\` object representing the volcano plot.
