% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readsPerMillion.R
\name{readsPerMillion}
\alias{readsPerMillion}
\title{Normalize Counts to Reads Per Million (RPM)}
\usage{
readsPerMillion(data, factor = 10^6)
}
\arguments{
\item{data}{A numeric matrix or data frame of raw counts. Rows are typically genes, and columns are samples.}

\item{factor}{A numeric value for scaling. Default is \code{10^6} to compute reads per million.}
}
\value{
A matrix of the same dimensions as \code{data}, with RPM-normalized values.
}
\description{
This function normalizes raw read counts to Reads Per Million (RPM) across each column (typically samples).
It is commonly used for comparing sequencing depth-normalized gene expression across samples.
}
