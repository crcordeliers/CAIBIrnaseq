#' Plot Expression Heatmap
#'
#' Generate a heatmap for gene expression data from a `SummarizedExperiment` object, with optional sample annotations and custom styling.
#'
#' @param expData A `SummarizedExperiment` object containing gene expression data.
#' @param genes A character vector of gene identifiers (e.g., `gene_name`) to include in the heatmap.
#' @param annotations A character vector of column names in the `colData` of `expData` to use for sample annotations in the heatmap. Default is `NA` (no annotations).
#' @param assay A character string specifying the assay to use from `expData`. Default is `"norm"`.
#' @param gene_name A character string specifying the gene identifier column in `rowData(expData)`. Default is `"gene_name"`.
#' @param annotation_prop A numeric value specifying the proportion of the heatmap height/width allocated to the annotation tracks. Default is `0.1`.
#' @param annotation_colors A named list of colors for the annotation tracks. Default is `NULL` (automatic coloring).
#' @param fname A character string specifying the file name for saving the heatmap. Default is `NULL` (no file saved).
#' @param fwidth Numeric, the width of the saved plot in inches. Default is `7`.
#' @param fheight Numeric, the height of the saved plot in inches. Default is `5`.
#' @param ... Additional arguments passed to the underlying heatmap plotting functions.
#'
#' @return A `ggplot` object representing the heatmap.
#'
#' @details
#' The function extracts and scales the expression data for the specified genes, optionally adds sample annotations from the `colData` of the `SummarizedExperiment`, and plots a heatmap using hierarchical clustering.
#' If a file name (`fname`) is provided, the heatmap is saved to the specified location.
#'
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom ggplot2 ggplot

#'
#' @export
#'
plot_exp_heatmap <- function(expData,
                             genes,
                             annotations = NA,
                             assay = "norm",
                             gene_name = "gene_name",
                             annotation_prop = 0.1,
                             annotation_colors = NULL,
                             fname = NULL,
                             fwidth = 7,
                             fheight = 5,
                             ...) {
  # Prepare the heatmap data (gene expression values)
  hm_data <- prep_exp_hm(expData, genes, assay, gene_name)

  # Generate the heatmap using the plt_heatmap function
  hm <- plt_heatmap(hm_data,
                    annotations = annotations,
                    fontsize = 8,
                    colors_title = "Scaled gene exp",
                    center = TRUE,
                    scale = TRUE,
                    hm_colors = 'RdBu',
                    raster = TRUE,
                    track_prop = annotation_prop,
                    track_colors = annotation_colors,
                    fname = fname,
                    fwidth = fwidth,
                    fheight = fheight,
                    ...)

  # Return the heatmap object
  return(hm)
}
