#' Save or return a scatter plot of two genes from a SummarizedExperiment
#'
#' This function extracts expression values for two genes from a `SummarizedExperiment` object,
#' generates a scatter plot, and optionally saves it to a file.
#'
#' @param exp_data A `SummarizedExperiment` object containing expression data and sample metadata.
#' @param gene1 A character string indicating the name of the gene for the x-axis.
#' @param gene2 A character string indicating the name of the gene for the y-axis.
#' @param color_var Optional character string indicating a column in metadata used for coloring points. Default is `NA`.
#' @param pt_size Numeric. Size of the points in the scatter plot. Default is 2.
#' @param fname Optional. A file name (with extension, e.g., `"scatterplot.png"`) to save the plot. If `NULL`, the plot is not saved.
#' @param fwidth Numeric. Width of the output figure in inches. Default is 5.
#' @param fheight Numeric. Height of the output figure in inches. Default is 3.
#'
#' @return A `ggplot2` object (the scatter plot).
#' @export
#'
#' @importFrom ggplot2 ggsave

plot_exp_scatter <- function(exp_data, gene1, gene2,
                             color_var = NA, pt_size = 2,
                             fname = NULL,
                             fwidth = 5,
                             fheight = 3) {
  exp_df <- get_exp_df(exp_data, c(gene1, gene2))
  plt <- plt_scatter(exp_df, gene1, gene2, color_var, pt_size)

  if (!is.null(fname)) {
    message("-- Saving plot at ", fname)
    ggplot2::ggsave(filename = fname, plot = plt, width = fwidth, height = fheight)
  }

  return(plt)
}
