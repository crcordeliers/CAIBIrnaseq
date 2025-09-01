#' Rebase Gene Expression Matrix
#'
#' This function rebases the gene expression matrix in a `SummarizedExperiment` object using a specified annotation as the primary identifier.
#'
#' @param exp_data A `SummarizedExperiment` object containing the gene expression
#' data.
#' @param gene_id A character string specifying the column in the gene annotation
#'  with the old main annotation. Default is `"gene_id"`.
#' @param new_gene_id A character string specifying the column in the gene annotation
#' to use as the main identifier. Default is `"gene_name"`.
#' @param keep_cols A character vector specifying columns to keep in the new
#' annotation. If `NULL`, all original columns are kept and aggregated as follows:
#' - For numeric variables, the mean is kept for multimapping resolution
#' - For character variables, all values are kept in a long string
#'
#' @return A `SummarizedExperiment` object with rebased gene expression data.
#' The output includes updated `rowData` with summarized gene-level information,
#' and new `assays` containing rebased `counts` and `tpm` matrices.
#'
#' @details
#' This function rebases the gene expression data by aggregating counts based on
#' the specified annotation. It also calculates transcripts per million (TPM)
#' using the aggregated data and associated gene lengths.
#'
#' The function performs the following steps:
#' - Aggregates gene counts and metadata based on the specified annotation.
#' - If `gene_lengths_kb` is available, calculates TPM values using the rebased
#' gene counts and average gene lengths.
#' - Constructs a new `SummarizedExperiment` object with the rebased data.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData colData assays SummarizedExperiment
#' @importFrom dplyr group_by summarize pull across
#' @importFrom rlang sym
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#'
rebase_gexp <- function(exp_data, gene_id = "gene_id", new_gene_id = "gene_name", keep_cols = NULL) {
  # Check input class
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("`exp_data` must be a SummarizedExperiment object.")
  }

  # Extract gene and sample annotations
  gene_annot <- as.data.frame(rowData(exp_data))
  sample_annot <- colData(exp_data)

  # Check if the annotation column exists
  if (!new_gene_id %in% colnames(gene_annot)) {
    stop("The specified annotation column '", new_gene_id, "' was not found in rowData.")
  }

  # Check for required columns
  required_cols <- c(gene_id, new_gene_id)
  missing_cols <- setdiff(required_cols, colnames(gene_annot))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s) in rowData: ", paste(missing_cols, collapse = ", "))
  }

  # Check that counts assay exists
  if (!"counts" %in% names(assays(exp_data))) {
    stop("Assay 'counts' not found in the SummarizedExperiment object.")
  }

  if(!all(rownames(assays(exp_data)[["counts"]]) %in% pull(gene_annot, !! gene_id))) {
    stop("Row names of the counts matrix must be found in gene_id: ",
         sum(rownames(assays(exp_data)[["counts"]]) %in% pull(gene_annot, !! gene_id)),
         " out of ",
         nrow(exp_data),
         " found.")
  }

  message("-- Rebasing the gene expression matrix using `", new_gene_id, "` as main annotation")

  # Preprocess gene expression matrix
  gexp_new <- gexp_preprocess(
    gexp = assays(exp_data)[["counts"]],
    gene_annotation = gene_annot,
    og_annot = gene_id,
    keep_annot = new_gene_id,
    keep_stat = "sum"
  )
  gexp_new <- as.matrix(gexp_new)


  if(is.null(keep_cols)) {
    keep_cols <- setdiff(colnames(gene_annot), c(gene_id, new_gene_id))
  }

  # Aggregate annotations
  new_gene_annot <- gene_annot |>
    group_by(!!sym(new_gene_id)) %>%
    summarize(across(.cols = where(is.numeric), .fns = mean),
              across(.cols = where(is.character), .fns = ~ paste0(unique(.), collapse = ", "))) |>
    select(!! new_gene_id, !! gene_id, !!! keep_cols)

  if("gene_length_kb" %in% colnames(new_gene_annot)) {
    message("Calculating new TPM...")
    lengths_kb <- pull(new_gene_annot, gene_length_kb, !!sym(new_gene_id))
    tpm <- transcripts_per_million(gexp_new, lengths_kb)

    new_exp_data <- SummarizedExperiment(
      assays = list(counts = gexp_new, tpm = tpm),
      colData = sample_annot,
      rowData = as_tibble(new_gene_annot)
    )
  } else {
    message("`gene_length_kb` not found in gene annotation: skipping TPM...")

    new_exp_data <- SummarizedExperiment(
      assays = list(counts = gexp_new),
      colData = sample_annot,
      rowData = as_tibble(new_gene_annot)
    )
  }

  return(new_exp_data)
}
