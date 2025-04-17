#' Extract gene expression and sample annotations
#'
#' This function extracts expression data for a given list of genes
#' from a SummarizedExperiment object and merges it with sample metadata.
#'
#' @param exp_data A `SummarizedExperiment` object.
#' @param genes A character vector of gene names or identifiers to extract.
#' @param assay A character string specifying which assay to extract. Defaults to `"norm"`.
#'
#' @return A data frame with expression values and sample annotations.
#' @export
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#'
get_exp_df <- function(exp_data, genes, assay = "norm") {
  gexp <- SummarizedExperiment::assays(exp_data)[[assay]]
  sample_annot <- SummarizedExperiment::colData(exp_data) |> as.data.frame()

  if (length(genes) > 1) {
    exp_df <- gexp[genes, ] |>
      t() |>
      as.data.frame()
  } else {
    exp_vec <- gexp[genes, ]
    exp_df <- matrix(
      exp_vec, ncol = 1,
      dimnames = list(colnames(gexp), genes)
    ) |> as.data.frame()
  }

  exp_df <- exp_df |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(sample_annot, by = "sample_id")

  return(exp_df)
}
