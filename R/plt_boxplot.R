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
#' @importFrom dplyr group_by summarize pull n
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin theme_light stat_summary
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
    group_by(!!sym(annotation)) |>
    summarize(count = n()) |>
    pull(count) |>
    max()

  use_box <- largest_n > 10 & summary_type != "line" | summary_type == "box"

  plt <- ggplot(exp_df, aes(!!sym(annotation), !!sym(gene)))

  # geom_violin needs at least 2 points per group to draw a density; below
  # that it silently drops the group, so only skip it for singleton groups.
  smallest_n <- exp_df |>
    group_by(!!sym(annotation)) |>
    summarize(count = n()) |>
    pull(count) |>
    min()

  if (smallest_n >= 2) {
    plt <- plt + geom_violin(aes(fill = !!sym(annotation)),
                              color = "grey40", alpha = 0.25,
                              trim = TRUE, show.legend = FALSE)
  }

  # Add beeswarm
  if (!is.na(color_var)) {
    plt <- plt + geom_beeswarm(aes(color = !!sym(color_var)), size = pt_size)
  } else {
    plt <- plt + geom_beeswarm(size = pt_size)
  }

  # Add summary (boxplot or mean line)
  if (use_box) {
    plt <- plt + geom_boxplot(width = 0.05, outlier.shape = NA)
  } else {
    plt <- plt + stat_summary(fun = mean, geom = "crossbar",
                                       fun.min = mean, fun.max = mean,
                                       width = 0.5, lwd = 0.2)
  }

  # Add stats
  if (is.na(stat_comparisons)) {
    plt <- plt + stat_compare_means(label = stat_format)
  } else {
    plt <- plt + stat_compare_means(comparisons = stat_comparisons, label = stat_format)
  }

  y_lab <- str_replace_all(gene, "_", " ")

  plt <- plt + theme_light() + labs(y = str_wrap(y_lab, width = 30))

  return(plt)
}
