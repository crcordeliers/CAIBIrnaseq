# Plot Heatmap with Optional Annotations and Customization

This function generates a heatmap from a data frame or matrix of values,
with the option to include additional annotations and customize various
aspects of the plot such as fonts, track colors, and size.

## Usage

``` r
plt_heatmap(
  data4hm,
  annotations = NA,
  fontsize = 9,
  track_colors = NULL,
  track_prop = 0.1,
  fname = NULL,
  fwidth = 7,
  fheight = 5,
  ...
)
```

## Arguments

- data4hm:

  A list containing the data for the heatmap with the several elements

- annotations:

  A data frame or list of annotations to overlay on the heatmap. Default
  is \`NA\`.

- fontsize:

  The font size for text in the heatmap. Default is 9.

- track_colors:

  A vector of colors for the annotation tracks. Default is \`NULL\`.

- track_prop:

  Proportion of the heatmap's height allocated to annotation tracks.
  Default is 0.1.

- fname:

  Optional. The file name to save the heatmap plot to. Default is
  \`NULL\`.

- fwidth:

  Width of the saved heatmap (in inches). Default is 7.

- fheight:

  Height of the saved heatmap (in inches). Default is 5.

- ...:

  Additional arguments to pass to the \`ggheatmap\` function.

## Value

A \`ggplot\` object representing the heatmap.
