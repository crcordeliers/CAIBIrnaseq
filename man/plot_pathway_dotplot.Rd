% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pathway_dotplot.R
\name{plot_pathway_dotplot}
\alias{plot_pathway_dotplot}
\title{Plot Pathway Analysis Results as a Dot Plot}
\usage{
plot_pathway_dotplot(
  exp_data,
  score_name = "resultsORA",
  top_n = 10,
  maxPval = 0.05
)
}
\arguments{
\item{exp_data}{A `SummarizedExperiment` object containing experimental data with pathway analysis results in `metadata()`.}

\item{score_name}{A character string indicating the metadata field where results are stored. Default is `"resultsORA"`.}

\item{top_n}{An integer for the number of top pathways to plot. Default is `10`.}

\item{maxPval}{A numeric for the maximum adjusted p-value to include pathways. Default is `0.05`.}
}
\value{
A `ggplot2` dot plot showing pathway enrichment results.
}
\description{
Creates a dot plot to visualize pathway analysis results, supporting both ORA and FGSEA formats.
}
\details{
The function supports two result formats:
- **ORA** (Over-Representation Analysis): expects columns `PAdj`, `GeneRatio`, and `Pathway`.
- **FGSEA**: expects columns `padj`, `size`, and `pathway`.

The x-axis shows the -log10 of adjusted p-values, y-axis lists pathways, dot size shows gene ratio.
}
