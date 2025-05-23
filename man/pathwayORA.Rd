% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathwayORA.R
\name{pathwayORA}
\alias{pathwayORA}
\title{Pathway Over-Representation Analysis (ORA)}
\usage{
pathwayORA(diffexp_result, pathways, id_col = "gene_symbol", pcutoff = 0.05)
}
\arguments{
\item{diffexp_result}{A data frame of differential expression results,
where row names correspond to gene identifiers.}

\item{pathways}{A data frame with at least two columns: one indicating
pathway names (e.g., 'pathway') and one with gene identifiers.}

\item{id_col}{Character. The column name in `pathways` that matches gene identifiers in `diffexp_result`.}

\item{pcutoff}{Numeric. Adjusted p-value threshold to filter significant pathways (default = 0.05).}
}
\value{
A data frame of enriched pathways with columns:
\describe{
  \item{Pathway}{Name of the enriched pathway}
  \item{PValue}{Raw p-value from Fisher's exact test}
  \item{PAdj}{Adjusted p-value (Benjamini-Hochberg)}
  \item{GeneRatio}{Proportion of input genes found in the pathway}
  \item{BgRatio}{Proportion of background genes found in the pathway}
  \item{Genes}{Comma-separated list of matched genes}
}
}
\description{
Performs an over-representation analysis (ORA) using Fisher's exact test
on a list of differentially expressed genes against a set of predefined pathways.
}
