#' Plot Expression Boxplot
#'
#' This function generates a boxplot for the expression of a specified gene, grouped by a given annotation.
#'
#' @importFrom ggplot2 ggsave
#' @param exp_data A data object containing expression data.
#' @param gene A character string specifying the gene for which expression is plotted.
#' @param annotation A character string specifying the grouping variable for the boxplot.
#' @param color_var Optional. A character string specifying the variable to color the points. Defaults to NA.
#' @param pt_size Numeric. The size of the points in the boxplot. Defaults to 1.
#' @param summary_type A character string specifying the type of summary to display ("choose", "line", or "box"). Defaults to "choose".
#' @param stat_comparisons Optional. Statistical comparisons to be displayed on the plot. Defaults to NA.
#' @param stat_format Optional. A format for statistical annotations. Defaults to NULL.
#' @param fname Optional. A character string specifying the file name to save the plot. Defaults to NULL.
#' @param fwidth Numeric. The width of the saved plot file. Defaults to 5.
#' @param fheight Numeric. The height of the saved plot file. Defaults to 3.
#'
#' @return A ggplot object representing the boxplot.
#'
#' @importFrom ggplot2 ggsave
#'
#' @export
#'
plot_exp_boxplot <- function(exp_data, gene, annotation,
                             color_var = NA,
                             pt_size = 1,
                             summary_type = c("choose", "line", "box")[1],
                             stat_comparisons = NA,
                             stat_format = NULL,
                             fname = NULL,
                             fwidth = 5,
                             fheight = 3) {
  exp_df <- get_exp_df(exp_data, gene)

  # Ensure annotation rownames match exp_df sample IDs
  if (!is.null(rownames(annotation)) && any(!exp_df$sample %in% rownames(annotation))) {
    stop("Some samples in expression data are missing in annotation data.")
  }

  plt <- plt_boxplot(
    exp_df = exp_df,
    gene = gene,
    annotation = annotation,
    color_var = color_var,
    pt_size = pt_size,
    summary_type = summary_type,
    stat_comparisons = stat_comparisons,
    stat_format = stat_format
  )

  if (!is.null(fname)) {
    message("-- Saving plot at ", fname)
    ggplot2::ggsave(
      filename = fname,
      plot = plt,
      width = fwidth,
      height = fheight,
      units = "in",
      dpi = 300
    )
  }

  return(plt)
}

