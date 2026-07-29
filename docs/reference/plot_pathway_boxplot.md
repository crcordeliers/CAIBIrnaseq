# Plot Pathway Boxplot

Creates a boxplot to visualize pathway scores across annotations, with
optional statistical comparisons and customization.

## Usage

``` r
plot_pathway_boxplot(
  exp_data,
  pathway,
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

  A \`SummarizedExperiment\` object containing pathway scores in the
  \`metadata(exp_data)\[\["pathway_scores"\]\]\` slot.

- pathway:

  A character string specifying the name of the pathway to plot.

- annotation:

  A character string specifying the column in \`colData(exp_data)\` to
  use for grouping samples in the plot.

- color_var:

  A character string specifying a column in \`colData(exp_data)\` to use
  for coloring the plot. Default is \`NA\` (no coloring).

- pt_size:

  Numeric; size of the points in the plot. Default is \`1\`.

- summary_type:

  Character; type of summary to overlay on the plot. Options are
  \`"choose"\` (default), \`"line"\`, or \`"box"\`.

- stat_comparisons:

  A list of character vectors specifying groups for pairwise statistical
  comparisons. Default is \`NA\` (no statistical comparisons).

- stat_format:

  A character string specifying the format for displaying statistical
  results. Default is \`NULL\` (no formatting).

- fname:

  A character string specifying the file name to save the plot. Default
  is \`NULL\` (do not save).

- fwidth:

  Numeric; width of the saved plot in inches. Default is \`5\`.

- fheight:

  Numeric; height of the saved plot in inches. Default is \`3\`.

## Value

A ggplot object representing the boxplot.

## Details

The function extracts pathway scores from the
\`metadata(exp_data)\[\["pathway_scores"\]\]\` slot and merges them with
sample annotations from \`colData(exp_data)\`. It generates a boxplot
using a helper function \`.plt_boxplot\`. Statistical comparisons can be
added to the plot, and the plot can be saved as an image file if
\`fname\` is provided.

The \`summary_type\` argument determines the type of summary overlay:

- "choose":

  No additional summary overlay.

- "line":

  Adds a line connecting the median of each group.

- "box":

  Adds a summary box around each group.
