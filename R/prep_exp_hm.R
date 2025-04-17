#' Prepare Gene Expression Data for Heatmap Visualization
#'
#' This function extracts and formats gene expression data from a `SummarizedExperiment`
#' object for heatmap plotting, returning a long-format data frame along with sample and feature axes.
#'
#' @param expData A `SummarizedExperiment` object.
#' @param genes Character vector of gene names or IDs to extract.
#' @param assay Character. The name of the assay to use (default is "norm").
#' @param gene_name Character. The column name in `rowData` containing gene names. Default is "gene_name".
#'
#' @returns A list containing:
#' \describe{
#'   \item{table}{A data frame with expression values and sample annotations.}
#'   \item{colv}{Column variable to use as x-axis (typically sample IDs).}
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
  gene_annot <- SummarizedExperiment::rowData(expData)

  if (!any(genes %in% rownames(expData)) && !any(genes %in% gene_annot[[gene_name]])) {
    stop("No `genes` were found in the row names or in the specified gene_name column of `expData`.")
  }

  # Match by gene_name column if needed
  if (!all(genes %in% rownames(expData)) && all(genes %in% gene_annot[[gene_name]])) {
    matched_genes <- gene_annot[[gene_name]] %in% genes
    rownames_to_use <- rownames(gene_annot)[matched_genes]
  } else {
    rownames_to_use <- genes[genes %in% rownames(expData)]
  }

  gexp <- SummarizedExperiment::assays(expData)[[assay]][rownames_to_use, , drop = FALSE]

  gexp_t <- gexp |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tibble::as_tibble()

  feats <- colnames(gexp_t)[-1]  # all but sample_id

  samp_annot <- SummarizedExperiment::colData(expData) |>
    as.data.frame()

  # Supprimer la colonne 'sample_id' si elle existe déjà
  if ("sample_id" %in% colnames(samp_annot)) samp_annot$sample_id <- NULL

  # Ajouter les rownames dans une colonne 'sample_id'
  samp_annot <- tibble::rownames_to_column(samp_annot, "sample_id")


  gexp_t <- dplyr::left_join(gexp_t, samp_annot, by = "sample_id")

  return(list(table = gexp_t, colv = "sample_id", rowv = feats))
}
