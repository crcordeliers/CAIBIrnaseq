# Boxplot of gene expression by annotation group

Creates a boxplot or crossbar summary of gene expression across sample
annotations, optionally colored and annotated with statistical
comparisons.

## Usage

``` r
plt_boxplot(
  exp_df,
  gene,
  annotation,
  color_var,
  pt_size,
  summary_type,
  stat_comparisons,
  stat_format
)
```

## Arguments

- exp_df:

  A data frame, typically created by \`get_exp_df()\`, containing gene
  expression and metadata.

- gene:

  A character string specifying the gene name to plot on the y-axis.

- annotation:

  A character string specifying the annotation variable (e.g.,
  "condition") on the x-axis.

- color_var:

  Optional character string specifying a variable for coloring the
  points. Set to \`NA\` to disable.

- pt_size:

  Numeric value for the size of the beeswarm points.

- summary_type:

  Character string: either \`"box"\` for boxplot or \`"line"\` for a
  mean crossbar.

- stat_comparisons:

  A list of comparison pairs for \`stat_compare_means()\`. Set to \`NA\`
  to use default behavior.

- stat_format:

  Character string for the label format used by \`stat_compare_means()\`
  (e.g., \`"p.signif"\`).

## Value

A \`ggplot2\` object representing the boxplot or summary plot.
