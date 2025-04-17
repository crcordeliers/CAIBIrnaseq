#' Boxplot of gene expression by annotation group
#'
#' Creates a boxplot or crossbar summary of gene expression across sample annotations,
#' optionally colored and annotated with statistical comparisons.
#'
#' @param exp_df A data frame, typically created by `get_exp_df()`, containing gene expression and metadata.
#' @param gene A character string specifying the gene name to plot on the y-axis.
#' @param annotation A character string specifying the annotation variable (e.g., "condition") on the x-axis.
#' @param color_var Optional character string specifying a variable for coloring the points. Set to `NA` to disable.
#' @param pt_size Numeric value for the size of the beeswarm points.
#' @param summary_type Character string: either `"box"` for boxplot or `"line"` for a mean crossbar.
#' @param stat_comparisons A list of comparison pairs for `stat_compare_means()`. Set to `NA` to use default behavior.
#' @param stat_format Character string for the label format used by `stat_compare_means()` (e.g., `"p.signif"`).
#'
#' @return A `ggplot2` object representing the boxplot or summary plot.
#' @export
#'
#' @importFrom dplyr group_by summarize pull
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_light stat_summary
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom ggpubr stat_compare_means
#' @importFrom rlang sym
#'
plt_boxplot <- function(exp_df, gene, annotation,
                        color_var,
                        pt_size,
                        summary_type,
                        stat_comparisons,
                        stat_format) {
  largest_n <- exp_df |>
    dplyr::group_by(!!rlang::sym(annotation)) |>
    dplyr::summarize(count = dplyr::n()) |>
    dplyr::pull(count) |>
    max()

  plt <- ggplot2::ggplot(exp_df, ggplot2::aes(!!rlang::sym(annotation), !!rlang::sym(gene)))

  # Add beeswarm
  if (!is.na(color_var)) {
    plt <- plt + ggbeeswarm::geom_beeswarm(ggplot2::aes(color = !!rlang::sym(color_var)), size = pt_size)
  } else {
    plt <- plt + ggbeeswarm::geom_beeswarm(size = pt_size)
  }

  # Add summary (boxplot or mean line)
  if (largest_n > 10 & summary_type != "line" | summary_type == "box") {
    plt <- plt + ggplot2::geom_boxplot(width = 0.05)
  } else {
    plt <- plt + ggplot2::stat_summary(fun = mean, geom = "crossbar",
                                       fun.min = mean, fun.max = mean,
                                       width = 0.5, lwd = 0.2)
  }

  # Add stats
  if (is.na(stat_comparisons)) {
    plt <- plt + ggpubr::stat_compare_means(label = stat_format)
  } else {
    plt <- plt + ggpubr::stat_compare_means(comparisons = stat_comparisons, label = stat_format)
  }

  plt <- plt + ggplot2::theme_light()
  return(plt)
}
