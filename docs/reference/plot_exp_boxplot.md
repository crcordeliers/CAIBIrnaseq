# Plot Expression Boxplot (with validation)

This function generates a boxplot for the expression of a specified
gene, grouped by a given annotation.

## Usage

``` r
plot_exp_boxplot(
  exp_data,
  gene,
  annotation,
  color_var = NA,
  pt_size = 1,
  summary_type = c("choose", "line", "box")[1],
  stat_comparisons = NA,
  stat_format = NULL,
  fname = NULL,
  fwidth = 5,
  fheight = 3
)
```

## Arguments

- exp_data:

  A data object containing expression data.

- gene:

  A character string specifying the gene for which expression is
  plotted.

- annotation:

  A character string specifying the grouping variable for the boxplot.

- color_var:

  Optional. A character string specifying the variable to color the
  points. Defaults to NA.

- pt_size:

  Numeric. The size of the points in the boxplot. Defaults to 1.

- summary_type:

  A character string specifying the type of summary to display
  ("choose", "line", or "box"). Defaults to "choose".

- stat_comparisons:

  Optional. Statistical comparisons to be displayed on the plot.
  Defaults to NA.

- stat_format:

  Optional. A format for statistical annotations. Defaults to NULL.

- fname:

  Optional. A character string specifying the file name to save the
  plot. Defaults to NULL.

- fwidth:

  Numeric. The width of the saved plot file. Defaults to 5.

- fheight:

  Numeric. The height of the saved plot file. Defaults to 3.

## Value

A ggplot object representing the boxplot.
