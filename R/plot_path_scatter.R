#' Plot a scatter plot of two pathway scores
#'
#' This function generates a scatter plot of two selected pathway scores extracted from
#' a `SummarizedExperiment` object. It uses the `pathway_scores` metadata and
#' optionally colors the points based on a variable from the sample annotations.
#'
#' @param exp_data A `SummarizedExperiment` object containing the pathway scores in its metadata.
#' @param pathway1 A character string indicating the name of the first pathway to plot.
#' @param pathway2 A character string indicating the name of the second pathway to plot.
#' @param color_var A character string indicating the column in the metadata for coloring the points. Default is `NA` (no coloring).
#' @param pt_size Numeric. The size of the points in the scatter plot. Default is 2.
#' @param fname Optional. A file name to save the plot (e.g., `"scatterplot.png"`). If `NULL`, the plot will not be saved.
#' @param fwidth Numeric. The width of the output figure in inches. Default is 5.
#' @param fheight Numeric. The height of the output figure in inches. Default is 3.
#'
#' @returns A `ggplot2` object containing the scatter plot.
#' @export
#'
#' @importFrom ggplot2 ggsave
#'
plot_path_scatter <- function(exp_data, pathway1, pathway2,
                              color_var = NA, pt_size = 2,
                              fname = NULL,
                              fwidth = 5,
                              fheight = 3){

  # Extract pathway scores from the metadata of the SummarizedExperiment object
  pathway_scores <- S4Vectors::metadata(exp_data)[["pathway_scores"]]

  # Check if pathway scores are present in the metadata
  if(is.null(pathway_scores)) {
    stop('No pathway scores found in `metadata(exp_data)[["pathway_scores"]]`')
  }

  # Get the pathway data frame for the selected pathways
  path_df <- get_pathway_df(exp_data, pathway_scores, c(pathway1, pathway2))

  # Create a scatter plot using the helper function `plt_scatter`
  plt <- plt_scatter(path_df, pathway1, pathway2, color_var, pt_size)

  # If a file name is provided, save the plot
  if(!is.null(fname)) {
    message("-- Saving plot at ", fname)
    ggplot2::ggsave(fname, plt, width = fwidth, height = fheight, create.dir = TRUE)
  }

  # Return the plot object
  return(plt)
}
