# Plot Pathway Scatter Plot

Creates a scatter plot comparing scores for two pathways from the
pathway scores stored in the metadata of a given dataset. Optionally
saves the plot to a file.

## Usage

``` r
plot_pathway_scatter(
  exp_data,
  pathway1,
  pathway2,
  color_var = NA,
  pt_size = 2,
  fname = NULL,
  fwidth = 5,
  fheight = 3
)
```

## Arguments

- exp_data:

  An object containing experimental data. Must include pathway scores in
  \`metadata(exp_data)\[\["pathway_scores"\]\]\`.

- pathway1:

  A string specifying the name of the first pathway for comparison.

- pathway2:

  A string specifying the name of the second pathway for comparison.

- color_var:

  (Optional) A variable used to color points in the scatter plot.
  Default is \`NA\`.

- pt_size:

  A numeric value specifying the size of points in the scatter plot.
  Default is \`2\`.

- fname:

  (Optional) A string specifying the file name to save the plot. If
  \`NULL\`, the plot is not saved. Default is \`NULL\`.

- fwidth:

  A numeric value specifying the width of the saved plot in inches.
  Default is \`5\`.

- fheight:

  A numeric value specifying the height of the saved plot in inches.
  Default is \`3\`.

## Value

A \`ggplot\` object representing the scatter plot.

## Details

This function extracts pathway scores from the metadata of the provided
experimental data object, then creates a scatter plot comparing the
scores of the specified pathways. An optional variable can be used to
color the points. If a file name is provided, the plot is saved to the
specified location.
