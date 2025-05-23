% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plt_scatter.R
\name{plt_scatter}
\alias{plt_scatter}
\title{Scatter plot of gene expression}
\usage{
plt_scatter(exp_df, gene1, gene2, color_var, pt_size)
}
\arguments{
\item{exp_df}{A data frame, typically created with `get_exp_df()`, containing expression values and metadata.}

\item{gene1}{A character string representing the first gene name (x-axis).}

\item{gene2}{A character string representing the second gene name (y-axis).}

\item{color_var}{(Optional) A character string naming the variable used for point coloring. Set to `NA` for no coloring.}

\item{pt_size}{Numeric. Size of the points in the plot.}
}
\value{
A `ggplot2` object representing the scatter plot with optional coloring and correlation statistics.
}
\description{
Creates a scatter plot of expression for two genes from a data frame, optionally colored by a variable.
}
