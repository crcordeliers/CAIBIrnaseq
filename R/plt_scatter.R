#' Scatter plot of gene expression
#'
#' Creates a scatter plot of expression for two genes from a data frame, optionally colored by a variable.
#'
#' @param exp_df A data frame, typically created with `get_exp_df()`, containing expression values and metadata.
#' @param gene1 A character string representing the first gene name (x-axis).
#' @param gene2 A character string representing the second gene name (y-axis).
#' @param color_var (Optional) A character string naming the variable used for point coloring. Set to `NA` for no coloring.
#' @param pt_size Numeric. Size of the points in the plot.
#'
#' @return A `ggplot2` object representing the scatter plot with optional coloring and correlation statistics.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_linedraw
#' @importFrom rlang sym
#' @importFrom ggpubr stat_cor
#'
plt_scatter <- function(exp_df, gene1, gene2, color_var, pt_size) {
  plt <- ggplot2::ggplot(exp_df, ggplot2::aes(!!rlang::sym(gene1), !!rlang::sym(gene2))) +
    ggpubr::stat_cor() +
    ggplot2::theme_linedraw()

  if (!is.na(color_var)) {
    plt <- plt + ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(color_var)), size = pt_size)
  } else {
    plt <- plt + ggplot2::geom_point(size = pt_size)
  }
  return(plt)
}
