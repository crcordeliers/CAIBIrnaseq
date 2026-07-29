# Plot Pathway Analysis Results as a Dot Plot

Creates a dot plot to visualize pathway enrichment analysis results,
supporting both ORA and FGSEA formats.

## Usage

``` r
plot_pathway_dotplot(exp_data, score_name = NULL, top_n = 10, maxPval = 0.05)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing pathway analysis results
  stored in the \`metadata()\`.

- score_name:

  A character string indicating the metadata field where pathway results
  are stored.

- top_n:

  Number of top pathways to display. Default is \`10\`.

- maxPval:

  Maximum adjusted p-value for filtering pathways (if \`usePval =
  TRUE\`). Default is \`0.05\`.

## Value

A \`ggplot2\` dot plot object showing pathway enrichment.

## Details

The function supports two types of enrichment results:

- ORA (Over-Representation Analysis): requires columns \`padj\`,
  \`geneRatio\`, and \`pathway\`

- FGSEA: requires columns \`padj\`, \`size\`, and \`pathway\`

When \`usePval = TRUE\`, the plot will show -log10 adjusted p-values on
the x-axis, colored by significance. When \`usePval = FALSE\`, the plot
will rank and size pathways by gene ratio or enrichment size.
