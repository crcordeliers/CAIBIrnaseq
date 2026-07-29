# Plot Microenvironment Heatmap

This function generates a heatmap of microenvironment population scores,
with optional sample annotations.

## Usage

``` r
plot_microenv_heatmap(
  exp_data,
  annotations = NA,
  annotation_prop = 0.1,
  annotation_colors = NULL,
  fname = NULL,
  fwidth = 7,
  fheight = 5,
  ...
)
```

## Arguments

- exp_data:

  A SummarizedExperiment object containing microenvironment population
  scores in its metadata.

- annotations:

  Optional. Sample annotations to include in the heatmap. Defaults to
  NA.

- annotation_prop:

  A numeric value specifying the proportion of the heatmap allocated to
  annotations. Defaults to 0.1.

- annotation_colors:

  Optional. A list of colors for annotations. Defaults to NULL.

- fname:

  Optional. A character string specifying the file name to save the
  heatmap. Defaults to NULL.

- fwidth:

  Numeric. The width of the saved heatmap file. Defaults to 7.

- fheight:

  Numeric. The height of the saved heatmap file. Defaults to 5.

- ...:

  Additional parameters passed to the heatmap plotting function.

## Value

A heatmap plot object.
