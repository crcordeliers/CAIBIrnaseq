# Plot Progeny Heatmap

Generates a heatmap for Progeny pathway scores, with optional sample
annotations and customization options.

## Usage

``` r
plot_progeny_heatmap(
  progeny_scores,
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

- progeny_scores:

  A data frame of pathway scores or a \`SummarizedExperiment\` object
  containing pathway scores in its metadata.

- annotations:

  (Optional) Annotations for the samples to be displayed as tracks on
  the heatmap. Default is \`NA\`, meaning no annotations are used.

- annotation_prop:

  Proportion of the heatmap width allocated to annotations, if provided.
  Default is \`0.1\`.

- annotation_colors:

  (Optional) A named list of colors for the annotation tracks. Names
  must match the annotation variables.

- fname:

  (Optional) File name to save the heatmap as an image. If \`NULL\`, the
  heatmap is displayed but not saved. Default is \`NULL\`.

- fwidth:

  Width of the saved heatmap image in inches. Default is \`7\`.

- fheight:

  Height of the saved heatmap image in inches. Default is \`5\`.

- ...:

  Additional arguments passed to the internal heatmap plotting function.

## Value

A heatmap object showing Progeny pathway scores.

## Details

If \`progeny_scores\` is a data frame, it is assumed to contain pathway
scores with row names as pathway names and column names as sample IDs.
If it is a \`SummarizedExperiment\` object, pathway scores are extracted
from the metadata under \`progeny_scores\`.

Annotations can be added to the heatmap if provided, and their
appearance can be customized using \`annotation_colors\`.
