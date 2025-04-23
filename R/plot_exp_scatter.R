#' Plot Expression Scatter Plot
#'
#' Creates a scatter plot to visualize the expression relationship between two genes, with optional coloring by a variable.
#'
#' @param exp_data A `SummarizedExperiment` object containing gene expression data in the assays slot.
#' @param gene1 A character string specifying the first gene for the x-axis.
#' @param gene2 A character string specifying the second gene for the y-axis.
#' @param color_var A character string specifying a column in `colData(exp_data)` to use for coloring the points. Default is `NA` (no coloring).
#' @param pt_size Numeric; size of the points in the scatter plot. Default is `2`.
#' @param fname A character string specifying the file name to save the plot. Default is `NULL` (do not save).
#' @param fwidth Numeric; width of the saved plot in inches. Default is `5`.
#' @param fheight Numeric; height of the saved plot in inches. Default is `3`.
#'
#' @details
#' The function extracts expression values for the specified genes (`gene1` and `gene2`) from the `assays(exp_data)` slot and combines them with sample annotations from `colData(exp_data)`.
#' It generates a scatter plot using a helper function `.plt_scatter`, with optional coloring based on a variable (`color_var`) from the sample annotations.
#'
#' The plot includes a correlation statistic computed with the `stat_cor()` function from the `ggpubr` package.
#'
#' If `fname` is provided, the plot will be saved as an image file.
#'
#' @return A ggplot object representing the scatter plot.
#'
#' @importFrom ggplot2 ggsave
#'
#' @export
#'
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
