#' Plot Expression Boxplot (with validation)
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
  # Basic validations
  if (missing(exp_data) || is.null(exp_data)) stop("Argument 'exp_data' is missing or NULL.")
  if (missing(gene) || is.null(gene) || !is.character(gene)) stop("Argument 'gene' must be a non-null character string.")
  if (missing(annotation) || is.null(annotation)) stop("Argument 'annotation' is missing or NULL.")

  if (!is.numeric(pt_size) || pt_size <= 0) stop("Argument 'pt_size' must be a positive number.")
  if (!summary_type %in% c("choose", "line", "box")) stop("Argument 'summary_type' must be one of 'choose', 'line', or 'box'.")
  if (!is.null(fname) && !is.character(fname)) stop("Argument 'fname' must be a character string if provided.")
  if (!is.numeric(fwidth) || fwidth <= 0) stop("Argument 'fwidth' must be a positive number.")
  if (!is.numeric(fheight) || fheight <= 0) stop("Argument 'fheight' must be a positive number.")

  exp_df <- get_exp_df(exp_data, gene)

  # Check if exp_df was successfully retrieved
  if (nrow(exp_df) == 0) {
    stop("No expression data retrieved for the specified gene: ", gene)
  }

  # Check if sample IDs match between expression and annotation
  if (!is.null(rownames(annotation)) && any(!exp_df$sample %in% rownames(annotation))) {
    stop("Some samples in expression data are missing in the annotation data.")
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
