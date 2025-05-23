% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_pathways.R
\name{score_pathways}
\alias{score_pathways}
\title{Score Pathways}
\usage{
score_pathways(exp_data, pathways, scoring_method = "gsva", verbose = TRUE)
}
\arguments{
\item{exp_data}{A `SummarizedExperiment` object containing normalized expression data in the `assays(exp_data)$norm` slot.}

\item{pathways}{A data frame with pathway definitions, containing at least two columns: `pathway` (pathway name) and either `gene_id` (Ensembl IDs) or `gene_symbol` (gene symbols).}

\item{scoring_method}{A character string specifying the scoring method to use. Options are `"gsva"`, `"ssgsea"`, `"plage"`, or `"zscore"`. Default is `"gsva"`.}

\item{verbose}{Logical; if `TRUE`, prints progress messages during computation. Default is `TRUE`.}
}
\value{
A data frame of pathway scores, with pathways as row names and samples as columns. Pathways are sorted by their total score variation.
}
\description{
Computes pathway scores for a given expression dataset using specified scoring methods.
}
\details{
The function identifies the gene annotation used in the expression matrix (`gene_id` or `gene_symbol`) by matching row names of `assays(exp_data)$norm` to the `pathways` data frame.
It then splits the pathways into gene sets and scores them using the specified method from the `GSVA` package.

The available scoring methods are:
\describe{
  \item{`gsva`}{Gene Set Variation Analysis.}
  \item{`ssgsea`}{Single-sample Gene Set Enrichment Analysis.}
  \item{`plage`}{Pathway Level Analysis of Gene Expression.}
  \item{`zscore`}{Z-score normalization.}
}

Pathway scores are sorted by the sum of absolute scores across samples, prioritizing pathways with the highest variation.
}
