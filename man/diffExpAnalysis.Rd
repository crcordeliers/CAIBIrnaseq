% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffExpAnalysis.R
\name{diffExpAnalysis}
\alias{diffExpAnalysis}
\title{Differential Expression Analysis}
\usage{
diffExpAnalysis(
  countData,
  sampleInfo,
  method = "DESeq2",
  cutoff = 10,
  annotation
)
}
\arguments{
\item{countData}{A matrix or data frame of raw count data. Rows represent genes, and columns represent samples.}

\item{sampleInfo}{A data frame containing sample metadata. Must include a `condition` column specifying the experimental conditions.}

\item{method}{A string specifying the method for differential expression analysis. Currently supports only `"DESeq2"`. Default is `"DESeq2"`.}

\item{cutoff}{An integer specifying the minimum number of counts required across all samples for a gene to be included in the analysis. Default is `10`.}

\item{annotation}{ee}
}
\value{
A data frame containing the results of the differential expression analysis, including adjusted p-values, log fold changes, and other statistics.
}
\description{
Performs differential expression analysis using DESeq2 on a given count matrix and sample information.
}
\details{
This function performs differential expression analysis using the DESeq2 package. It filters genes with low counts, estimates size factors for normalization, and performs the DESeq2 analysis pipeline. Log fold-change shrinkage is applied using the `lfcShrink` function.
}
