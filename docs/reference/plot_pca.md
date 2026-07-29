# Plot PCA Results for Gene Expression Data

This function generates a PCA plot for a \`SummarizedExperiment\` object
using precomputed PCA results. It creates a scatter plot of the
specified principal components (\`pcs\[1\]\` and \`pcs\[2\]\`). Metadata
can be visualized via point color and shape.

## Usage

``` r
plot_pca(
  exp_data,
  color = NA,
  shape = NA,
  fname = "results/qc/plot_PCA.pdf",
  pcs = c(1, 2),
  res_name = "pca_res",
  id_name = "sample_id",
  ellipses = FALSE,
  point_size = 2,
  out = c("plotly", "ggplot")[1]
)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing gene expression data. PCA
  results must be stored in \`exp_data@metadata\` under the name
  specified by \`res_name\`.

- color:

  A character string specifying the column name in \`colData(exp_data)\`
  to use for coloring the points. Default is \`NA\`, meaning no
  coloring.

- shape:

  A character string specifying the column name in \`colData(exp_data)\`
  to use for shaping the points. Default is \`NA\`, meaning no shaping.

- fname:

  A character string specifying the file path to save the plot as a PDF.
  Default is \`"results/qc/plot_PCA.pdf"\`. Set to \`NULL\` to skip
  saving.

- pcs:

  An integer vector of length 2 specifying the principal components to
  plot on the x and y axes. Default is \`c(1, 2)\`.

- res_name:

  A character string specifying the name of the PCA results stored in
  \`exp_data@metadata\`. Default is \`"pca_res"\`.

- id_name:

  A character string specifying the column name in \`colData(exp_data)\`
  containing sample identifiers. Default is \`"sample_id"\`.

- ellipses:

  if TRUE, \`stat_ellipse\` is added to the points, representing the
  0.95 confidence-interval ellipse for each group in \`color\`.

- point_size:

  A numeric value specifying the size of the points in the plot. Default
  is \`2\`.

- out:

  A character string indicating the output type: \`"plotly"\`
  (interactive Plotly plot) or \`"ggplot"\` (static ggplot). Default is
  \`"plotly"\`.

## Value

A PCA plot, either as a \`plotly\` interactive object or a \`ggplot\`
static object, depending on the \`out\` parameter.
