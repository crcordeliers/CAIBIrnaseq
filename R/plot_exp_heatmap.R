#' Plot Heatmap of Gene Expression with Optional Annotations
#'
#' This function generates a heatmap of gene expression data from a `SummarizedExperiment` object,
#' with the option to add annotations to the heatmap. The data is pre-processed before plotting,
#' and the heatmap can be customized in terms of colors, font size, and file output.
#'
#' @param expData A `SummarizedExperiment` object containing the expression data. The `norm` assay is used by default.
#' @param genes A character vector of gene names or IDs to include in the heatmap.
#' @param annotations A data frame or list of annotations to overlay on the heatmap. Default is `NA` (no annotations).
#' @param assay The assay to use from the `expData` object. Default is `"norm"`.
#' @param gene_name The column name in `rowData(expData)` that contains the gene identifiers. Default is `"gene_name"`.
#' @param annotation_prop The proportion of the heatmap's height to allocate to the annotation tracks. Default is 0.1.
#' @param annotation_colors A vector of colors for the annotation tracks. Default is `NULL`.
#' @param fname Optional. The file name to save the heatmap plot to. Default is `NULL` (no file is saved).
#' @param fwidth The width (in inches) of the saved heatmap. Default is 7 inches.
#' @param fheight The height (in inches) of the saved heatmap. Default is 5 inches.
#' @param ... Additional arguments to pass to the `plt_heatmap` function.
#'
#' @returns A `ggplot` object representing the heatmap.
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
