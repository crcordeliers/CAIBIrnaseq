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
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom ggplot2 ggplot
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

  if (!inherits(expData, "SummarizedExperiment")) {
    stop("expData must be a SummarizedExperiment object.")
  }
  if (missing(genes) || !is.character(genes) || length(genes) == 0) {
    stop("You must provide a non-empty character vector for genes.")
  }
  if (!is.character(assay) || length(assay) != 1) {
    stop("assay must be a single character string.")
  }
  if (!assay %in% SummarizedExperiment::assayNames(expData)) {
    stop(paste0("Assay '", assay, "' not found in expData. Available assays are: ",
                paste(SummarizedExperiment::assayNames(expData), collapse = ", ")))
  }
  if (!is.character(gene_name) || length(gene_name) != 1) {
    stop("gene_name must be a single character string.")
  }
  if (!gene_name %in% colnames(SummarizedExperiment::rowData(expData))) {
    stop(paste0("gene_name '", gene_name, "' not found in rowData(expData)."))
  }
  if (!is.numeric(annotation_prop) || annotation_prop < 0 || annotation_prop > 1) {
    stop("annotation_prop must be a numeric value between 0 and 1.")
  }
  if (!is.null(fname) && !is.character(fname)) {
    stop("fname must be NULL or a character string (path to save the figure).")
  }

  ## Optional: validate annotations if provided
  if (!is.na(annotations[1])) {
    if (!all(annotations %in% colnames(SummarizedExperiment::colData(expData)))) {
      missing_annots <- annotations[!annotations %in% colnames(SummarizedExperiment::colData(expData))]
      stop(paste0("Some annotations not found in expData: ", paste(missing_annots, collapse = ", ")))
    }
  }

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

  return(hm)
}
