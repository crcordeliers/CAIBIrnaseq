#' Prepare Gene Expression Data for Heatmap Visualization
#'
#' This function extracts and formats gene expression data from a `SummarizedExperiment`
#' object for heatmap plotting. It returns a long-format data frame along with identifiers for samples and features.
#'
#' @param expData A `SummarizedExperiment` object.
#' @param genes Character vector of gene names or IDs to extract.
#' @param assay Character. The name of the assay to use (default is "norm").
#' @param gene_name Character. The column name in `rowData` containing gene names. Default is "gene_name".
#'
#' @returns A list containing:
#' \describe{
#'   \item{table}{A data frame with expression values and sample annotations.}
#'   \item{colv}{Column variable to use as the x-axis (typically sample IDs).}
#'   \item{rowv}{Row variable(s) corresponding to the genes.}
#' }
#' @export
#'
#' @importFrom SummarizedExperiment rowData assays colData
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr left_join
#'
prep_exp_hm <- function(expData,
                        genes,
                        assay = "norm",
                        gene_name = "gene_name") {
  # Basic input checks
  if (!inherits(expData, "SummarizedExperiment")) {
    stop("`expData` must be a 'SummarizedExperiment' object.")
  }

  if (!is.character(genes) || length(genes) == 0) {
    stop("`genes` must be a non-empty character vector.")
  }

  if (!is.character(assay) || length(assay) != 1) {
    stop("`assay` must be a character string of length 1.")
  }

  if (!(assay %in% names(SummarizedExperiment::assays(expData)))) {
    stop(paste0("Assay '", assay, "' not found in `expData`."))
  }

  if (!is.character(gene_name) || length(gene_name) != 1) {
    stop("`gene_name` must be a character string of length 1.")
  }

  gene_annot <- SummarizedExperiment::rowData(expData)

  if (!(gene_name %in% colnames(gene_annot))) {
    stop(paste0("Column `gene_name` '", gene_name, "' not found in `rowData(expData)`."))
  }

  if (!any(genes %in% rownames(expData)) && !any(genes %in% gene_annot[[gene_name]])) {
    stop("None of the specified `genes` were found in row names or the `gene_name` column of `expData`.")
  }

  # Determine genes to extract
  if (!all(genes %in% rownames(expData)) && all(genes %in% gene_annot[[gene_name]])) {
    matched_genes <- gene_annot[[gene_name]] %in% genes
    rownames_to_use <- rownames(gene_annot)[matched_genes]
  } else {
    rownames_to_use <- genes[genes %in% rownames(expData)]
  }

  if (length(rownames_to_use) == 0) {
    stop("No valid genes to extract after filtering.")
  }

  gexp <- SummarizedExperiment::assays(expData)[[assay]][rownames_to_use, , drop = FALSE]

  if (!is.matrix(gexp)) {
    stop("Selected assay must be a matrix.")
  }

  # Transform and join with sample annotations
  gexp_t <- gexp |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tibble::as_tibble()

  feats <- colnames(gexp_t)[-1]  # All but sample_id

  samp_annot <- SummarizedExperiment::colData(expData) |>
    as.data.frame()

  if ("sample_id" %in% colnames(samp_annot)) samp_annot$sample_id <- NULL

  samp_annot <- tibble::rownames_to_column(samp_annot, "sample_id")

  gexp_t <- dplyr::left_join(gexp_t, samp_annot, by = "sample_id")

  return(list(table = gexp_t, colv = "sample_id", rowv = feats))
}
