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
#' @importFrom tibble rownames_to_column as_tibble tibble
#' @importFrom dplyr left_join
#'
prep_exp_hm <- function(expData,
                        genes,
                        assay = "norm",
                        gene_name = "gene_name") {
  gene_annot <- SummarizedExperiment::rowData(expData)
  if (!any(genes %in% rownames(expData))) {
    stop("No `genes` were found as `gene_name` in the `exp_data` object.")
  } else if(!all(genes %in% gene_annot$gene_name)) {
    genes <- genes[genes %in% rownames(expData)]
  } else {
    genes <- genes
  }
  gexp <- SummarizedExperiment::assays(expData)[[assay]][genes,]
  gexp_t <- gexp |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tibble::tibble()

  feats <- colnames(gexp_t)[-1]

  samp_annot <- SummarizedExperiment::colData(exp_data) |> as.data.frame()

  gexp_t <- gexp_t |>
    dplyr::left_join(samp_annot, by = "sample_id")

  return(list(table = gexp_t, colv = "sample_id", rowv = feats))
}

