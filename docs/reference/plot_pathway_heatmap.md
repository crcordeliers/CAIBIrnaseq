# Plot Pathway Heatmap

This function generates a heatmap of pathway scores, with optional
sample annotations.

## Usage

``` r
plot_pathway_heatmap(
  pathway_scores,
  pathways = 1:20,
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

- pathway_scores:

  Either a data frame of pathway scores or a SummarizedExperiment object
  containing pathway scores in its metadata.

- pathways:

  A vector specifying the pathways to include in the heatmap. Defaults
  to the first 20 pathways.

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
