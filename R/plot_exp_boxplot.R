#' Plot expression boxplot for a gene
#'
#' This function generates a boxplot (or similar summary plot) of gene expression
#' using a provided expression matrix, gene of interest, and sample annotations.
#'
#' @param exp_data A data frame or matrix containing expression data with genes as rows and samples as columns.
#' @param gene Character. The gene name to be plotted.
#' @param annotation A data frame with sample annotations. Must contain row names or a column matching the sample names in `exp_data`.
#' @param color_var Character. Column name in `annotation` to use for color grouping. Default is `NA` (no coloring).
#' @param pt_size Numeric. Size of the individual points on the plot. Default is 1.
#' @param summary_type Character. Type of summary to add: "choose" (automatic), "line", or "box". Default is "choose".
#' @param stat_comparisons A list of group comparisons for statistical testing (passed to `ggpubr::stat_compare_means`). Default is `NA`.
#' @param stat_format A function or format for customizing statistical annotations. Default is `NULL`.
#' @param fname Character. If not `NULL`, the filename to save the plot to (PDF, PNG, etc.). Default is `NULL`.
#' @param fwidth Numeric. Width of the saved figure in inches. Default is 5.
#' @param fheight Numeric. Height of the saved figure in inches. Default is 3.
#'
#' @return A ggplot object containing the boxplot.
#' @export
#'
#' @importFrom ggplot2 ggsave
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

