#' Plot a boxplot of pathway scores by annotation group
#'
#' This function generates a boxplot (or alternative summaries) of pathway scores for a given pathway,
#' grouped by a chosen annotation variable (e.g., condition). It can optionally color the points,
#' and save the figure.
#'
#' @param exp_data A `SummarizedExperiment` object containing the pathway scores in its metadata.
#' @param pathway A character string indicating the name of the pathway to plot.
#' @param annotation A character string indicating the column in the sample metadata to group the data by (x-axis).
#' @param color_var A character string for the column to use for coloring points. Default is `NA`.
#' @param pt_size Numeric. Point size for the scatter/boxplot. Default is 1.
#' @param summary_type Character. One of `"box"`, `"line"`, or `"choose"` (default).
#' @param stat_comparisons Optional. Statistical comparison groups (for adding p-values). Default is `NA`.
#' @param stat_format Optional. Format string for displaying p-values. Default is `NULL`.
#' @param fname Optional. File name to save the figure. Default is `NULL` (no saving).
#' @param fwidth Numeric. Width of the saved plot in inches. Default is 5.
#' @param fheight Numeric. Height of the saved plot in inches. Default is 3.
#'
#' @return A `ggplot2` object representing the plot.
#' @export
#'
#' @importFrom ggplot2 ggsave
#'
plot_path_boxplot <- function(exp_data, pathway, annotation,
                              color_var = NA,
                              pt_size = 1,
                              summary_type = c("choose", "line", "box")[1],
                              stat_comparisons = NA,
                              stat_format = NULL,
                              fname = NULL,
                              fwidth = 5,
                              fheight = 3) {
  pathway_scores <- metadata(exp_data)[["pathway_scores"]]

  if (is.null(pathway_scores)) {
    stop('No pathway scores found in `metadata(exp_data)[["pathway_scores"]]`')
  }

  path_df <- get_pathway_df(exp_data, pathway_scores, pathway)

  plt <- plt_boxplot(path_df, pathway, annotation, color_var,
                     pt_size, summary_type, stat_comparisons, stat_format)

  if (!is.null(fname)) {
    message("-- Saving plot at ", fname)
    ggplot2::ggsave(fname, plt, width = fwidth, height = fheight)
  }

  return(plt)
}
