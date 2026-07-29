# Plot Expression Heatmap

Generate a heatmap for gene expression data from a
\`SummarizedExperiment\` object, with optional sample annotations and
custom styling.

## Usage

``` r
plot_exp_heatmap(
  expData,
  genes,
  annotations = NA,
  assay = "norm",
  gene_name = "gene_name",
  annotation_prop = 0.1,
  annotation_colors = NULL,
  fname = NULL,
  fwidth = 7,
  fheight = 5,
  ...
)
```

## Arguments

- expData:

  A \`SummarizedExperiment\` object containing gene expression data.

- genes:

  A character vector of gene identifiers (e.g., \`gene_name\`) to
  include in the heatmap.

- annotations:

  A character vector of column names in the \`colData\` of \`expData\`
  to use for sample annotations in the heatmap. Default is \`NA\` (no
  annotations).

- assay:

  A character string specifying the assay to use from \`expData\`.
  Default is \`"norm"\`.

- gene_name:

  A character string specifying the gene identifier column in
  \`rowData(expData)\`. Default is \`"gene_name"\`.

- annotation_prop:

  A numeric value specifying the proportion of the heatmap height/width
  allocated to the annotation tracks. Default is \`0.1\`.

- annotation_colors:

  A named list of colors for the annotation tracks. Default is \`NULL\`
  (automatic coloring).

- fname:

  A character string specifying the file name for saving the heatmap.
  Default is \`NULL\` (no file saved).

- fwidth:

  Numeric, the width of the saved plot in inches. Default is \`7\`.

- fheight:

  Numeric, the height of the saved plot in inches. Default is \`5\`.

- ...:

  Additional arguments passed to the underlying heatmap plotting
  functions.

## Value

A \`ggplot\` object representing the heatmap.
