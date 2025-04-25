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
get_exp_df <- function(exp_data, genes, assay = "norm") {
  # --- Input checks ---
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("`exp_data` must be a SummarizedExperiment object.")
  }

  if (!is.character(genes) || length(genes) == 0) {
    stop("`genes` must be a non-empty character vector.")
  }

  if (!assay %in% names(SummarizedExperiment::assays(exp_data))) {
    stop(paste0("Assay '", assay, "' not found in SummarizedExperiment object. Available assays: ",
                paste(names(SummarizedExperiment::assays(exp_data)), collapse = ", ")))
  }

  gexp <- SummarizedExperiment::assays(exp_data)[[assay]]
  sample_annot <- SummarizedExperiment::colData(exp_data) |> as.data.frame()

  # --- Validate gene presence ---
  missing_genes <- setdiff(genes, rownames(gexp))
  if (length(missing_genes) > 0) {
    warning("The following genes are not found in the assay and will be ignored: ",
            paste(missing_genes, collapse = ", "))
    genes <- intersect(genes, rownames(gexp))
  }

  if (length(genes) == 0) {
    stop("None of the specified genes were found in the assay.")
  }

  # --- Extract expression matrix ---
  if (length(genes) > 1) {
    exp_df <- gexp[genes, , drop = FALSE] |>
      t() |>
      as.data.frame()
  } else {
    exp_vec <- gexp[genes, ]
    exp_df <- matrix(
      exp_vec, ncol = 1,
      dimnames = list(colnames(gexp), genes)
    ) |> as.data.frame()
  }

  # --- Merge with sample metadata ---
  exp_df <- exp_df |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(sample_annot, by = "sample_id")

  return(exp_df)
}
