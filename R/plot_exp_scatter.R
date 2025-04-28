#' Plot Expression Scatter Plot
#'
#' Creates a scatter plot to visualize the expression relationship between two genes, with optional coloring by a variable.
#'
#' @param exp_data A `SummarizedExperiment` object containing gene expression data in the assays slot.
#' @param gene1 A character string specifying the first gene for the x-axis.
#' @param gene2 A character string specifying the second gene for the y-axis.
#' @param color_var A character string specifying a column in `colData(exp_data)` to use for coloring the points. Default is `NA` (no coloring).
#' @param pt_size Numeric; size of the points in the scatter plot. Default is `2`.
#' @param fname A character string specifying the file name to save the plot. Default is `NULL` (do not save).
#' @param fwidth Numeric; width of the saved plot in inches. Default is `5`.
#' @param fheight Numeric; height of the saved plot in inches. Default is `3`.
#'
#' @return A ggplot object representing the scatter plot.
#'
#' @importFrom ggplot2 ggsave
#' @export
#'
plot_exp_scatter <- function(exp_data, gene1, gene2,
                             color_var = NA, pt_size = 2,
                             fname = NULL,
                             fwidth = 5,
                             fheight = 3) {
  ## --------------------------- ##
  ## Verifications / Validation  ##
  ## --------------------------- ##
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("exp_data must be a SummarizedExperiment object.")
  }
  if (missing(gene1) || !is.character(gene1) || length(gene1) != 1) {
    stop("gene1 must be a single character string.")
  }
  if (missing(gene2) || !is.character(gene2) || length(gene2) != 1) {
    stop("gene2 must be a single character string.")
  }

  # Check if gene1 and gene2 exist in rowData
  gene_names <- SummarizedExperiment::rowData(exp_data)[["gene_name"]]
  if (is.null(gene_names)) {
    stop("No 'gene_name' column found in rowData(exp_data).")
  }
  if (!(gene1 %in% gene_names)) {
    stop(paste0("Gene '", gene1, "' not found in rowData(exp_data)."))
  }
  if (!(gene2 %in% gene_names)) {
    stop(paste0("Gene '", gene2, "' not found in rowData(exp_data)."))
  }

  if (!is.na(color_var)) {
    if (!is.character(color_var) || !(color_var %in% colnames(SummarizedExperiment::colData(exp_data)))) {
      stop(paste0("color_var '", color_var, "' not found in colData(exp_data)."))
    }
  }
  if (!is.numeric(pt_size) || pt_size <= 0) {
    stop("pt_size must be a positive numeric value.")
  }
  if (!is.null(fname) && !is.character(fname)) {
    stop("fname must be NULL or a character string (path to save the figure).")
  }

  ## --------------------------- ##
  ##         Main Function        ##
  ## --------------------------- ##

  exp_df <- get_exp_df(exp_data, c(gene1, gene2))
  plt <- plt_scatter(exp_df, gene1, gene2, color_var, pt_size)

  if (!is.null(fname)) {
    message("-- Saving plot at ", fname)
    ggplot2::ggsave(filename = fname, plot = plt, width = fwidth, height = fheight, dpi = 300)
  }

  return(plt)
}
