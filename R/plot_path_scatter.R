#' Plot Pathway Scatter Plot
#'
#' Creates a scatter plot comparing scores for two pathways from the pathway scores
#' stored in the metadata of a given dataset. Optionally saves the plot to a file.
#'
#' @param exp_data An object containing experimental data. Must include pathway scores in `metadata(exp_data)[["pathway_scores"]]`.
#' @param pathway1 A string specifying the name of the first pathway for comparison.
#' @param pathway2 A string specifying the name of the second pathway for comparison.
#' @param color_var (Optional) A variable used to color points in the scatter plot. Default is `NA`.
#' @param pt_size A numeric value specifying the size of points in the scatter plot. Default is `2`.
#' @param fname (Optional) A string specifying the file name to save the plot. If `NULL`, the plot is not saved. Default is `NULL`.
#' @param fwidth A numeric value specifying the width of the saved plot in inches. Default is `5`.
#' @param fheight A numeric value specifying the height of the saved plot in inches. Default is `3`.
#'
#' @details
#' This function extracts pathway scores from the metadata of the provided experimental
#' data object, then creates a scatter plot comparing the scores of the specified pathways.
#' An optional variable can be used to color the points. If a file name is provided, the
#' plot is saved to the specified location.
#'
#' @return A `ggplot` object representing the scatter plot.
#'
#' @importFrom S4Vectors metadata
#' @importFrom ggplot2 ggsave
#'
#' @export
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
