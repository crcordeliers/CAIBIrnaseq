# Scatter plot of gene expression

Creates a scatter plot of expression for two genes from a data frame,
optionally colored by a variable.

## Usage

``` r
plt_scatter(exp_df, gene1, gene2, color_var, pt_size)
```

## Arguments

- exp_df:

  A data frame, typically created with \`get_exp_df()\`, containing
  expression values and metadata.

- gene1:

  A character string representing the first gene name (x-axis).

- gene2:

  A character string representing the second gene name (y-axis).

- color_var:

  (Optional) A character string naming the variable used for point
  coloring. Set to \`NA\` for no coloring.

- pt_size:

  Numeric. Size of the points in the plot.

## Value

A \`ggplot2\` object representing the scatter plot with optional
coloring and correlation statistics.
