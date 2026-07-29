# Plot Expression Scatter Plot

Creates a scatter plot to visualize the expression relationship between
two genes, with optional coloring by a variable.

## Usage

``` r
plot_exp_scatter(
  exp_data,
  gene1,
  gene2,
  color_var = NA,
  pt_size = 2,
  fname = NULL,
  fwidth = 5,
  fheight = 3
)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing gene expression data in
  the assays slot.

- gene1:

  A character string specifying the first gene for the x-axis.

- gene2:

  A character string specifying the second gene for the y-axis.

- color_var:

  A character string specifying a column in \`colData(exp_data)\` to use
  for coloring the points. Default is \`NA\` (no coloring).

- pt_size:

  Numeric; size of the points in the scatter plot. Default is \`2\`.

- fname:

  A character string specifying the file name to save the plot. Default
  is \`NULL\` (do not save).

- fwidth:

  Numeric; width of the saved plot in inches. Default is \`5\`.

- fheight:

  Numeric; height of the saved plot in inches. Default is \`3\`.

## Value

A ggplot object representing the scatter plot.
