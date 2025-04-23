#' Plot Pathway Boxplot
#'
#' Creates a boxplot to visualize pathway scores across annotations, with optional statistical comparisons and customization.
#'
#' @param exp_data A `SummarizedExperiment` object containing pathway scores in the `metadata(exp_data)[["pathway_scores"]]` slot.
#' @param pathway A character string specifying the name of the pathway to plot.
#' @param annotation A character string specifying the column in `colData(exp_data)` to use for grouping samples in the plot.
#' @param color_var A character string specifying a column in `colData(exp_data)` to use for coloring the plot. Default is `NA` (no coloring).
#' @param pt_size Numeric; size of the points in the plot. Default is `1`.
#' @param summary_type Character; type of summary to overlay on the plot. Options are `"choose"` (default), `"line"`, or `"box"`.
#' @param stat_comparisons A list of character vectors specifying groups for pairwise statistical comparisons. Default is `NA` (no statistical comparisons).
#' @param stat_format A character string specifying the format for displaying statistical results. Default is `NULL` (no formatting).
#' @param fname A character string specifying the file name to save the plot. Default is `NULL` (do not save).
#' @param fwidth Numeric; width of the saved plot in inches. Default is `5`.
#' @param fheight Numeric; height of the saved plot in inches. Default is `3`.
#'
#' @details
#' The function extracts pathway scores from the `metadata(exp_data)[["pathway_scores"]]` slot and merges them with sample annotations from `colData(exp_data)`.
#' It generates a boxplot using a helper function `.plt_boxplot`. Statistical comparisons can be added to the plot, and the plot can be saved as an image file if `fname` is provided.
#'
#' The `summary_type` argument determines the type of summary overlay:
#' \describe{
#'   \item{"choose"}{No additional summary overlay.}
#'   \item{"line"}{Adds a line connecting the median of each group.}
#'   \item{"box"}{Adds a summary box around each group.}
#' }
#'
#' @return A ggplot object representing the boxplot.
#'
#' @importFrom S4Vectors metadata
#' @importFrom ggplot2 ggsave
#'
#' @export
plot_path_boxplot <- function(exp_data, pathway, annotation,
                              color_var = NA,
                              pt_size = 1,
                              summary_type = c("choose", "line", "box")[1],
                              stat_comparisons = NA,
                              stat_format = NULL,
                              fname = NULL,
                              fwidth = 5,
                              fheight = 3) {
  pathway_scores <- S4Vectors::metadata(exp_data)[["pathway_scores"]]

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
