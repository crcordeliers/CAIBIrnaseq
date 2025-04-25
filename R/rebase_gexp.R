#' Rebase Gene Expression Matrix
#'
#' This function rebases the gene expression matrix in a `SummarizedExperiment` object using a specified annotation as the primary identifier.
#'
#' @param exp_data A `SummarizedExperiment` object containing the gene expression data.
#' @param annotation A character string specifying the column in the gene annotation to use as the main identifier. Default is `"gene_name"`.
#'
#' @return A `SummarizedExperiment` object with rebased gene expression data. The output includes updated `rowData` with summarized gene-level information, and new `assays` containing rebased `counts` and `tpm` matrices.
#'
#' @details
#' This function rebases the gene expression data by aggregating counts based on the specified annotation. It also calculates transcripts per million (TPM) using the aggregated data and associated gene lengths.
#'
#' The function performs the following steps:
#' - Aggregates gene counts and metadata based on the specified annotation.
#' - Calculates TPM values using the rebased gene counts and average gene lengths.
#' - Constructs a new `SummarizedExperiment` object with the rebased data.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData colData assays SummarizedExperiment
#' @importFrom dplyr group_by summarize pull
#' @importFrom rlang sym
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#'
rebase_gexp <- function(exp_data, annotation = "gene_name") {
  # Check input class
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("`exp_data` must be a SummarizedExperiment object.")
  }

  # Extract gene and sample annotations
  gene_annot <- as.data.frame(SummarizedExperiment::rowData(exp_data))
  sample_annot <- SummarizedExperiment::colData(exp_data)

  # Check if the annotation column exists
  if (!annotation %in% colnames(gene_annot)) {
    stop("The specified annotation column '", annotation, "' was not found in rowData.")
  }

  # Check for required columns
  required_cols <- c("gene_id", "gene_length_kb", "gene_description", "gene_biotype")
  missing_cols <- setdiff(required_cols, colnames(gene_annot))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s) in rowData: ", paste(missing_cols, collapse = ", "))
  }

  # Check that counts assay exists
  if (!"counts" %in% names(SummarizedExperiment::assays(exp_data))) {
    stop("Assay 'counts' not found in the SummarizedExperiment object.")
  }

  message("-- Rebasing the gene expression matrix using `", annotation, "` as main annotation")

  # Preprocess gene expression matrix
  gexp_new <- gexp_preprocess(
    gexp = SummarizedExperiment::assays(exp_data)[["counts"]],
    gene_annotation = gene_annot,
    og_annot = "gene_id",
    keep_annot = annotation,
    keep_stat = "sum"
  )
  gexp_new <- as.matrix(gexp_new)

  # Aggregate annotations
  new_gene_annot <- gene_annot %>%
    dplyr::group_by(!!sym(annotation)) %>%
    dplyr::summarize(
      gene_id = paste(unique(gene_id), collapse = ", "),
      gene_length_kb = mean(gene_length_kb, na.rm = TRUE),
      gene_description = paste(unique(gene_description), collapse = ", "),
      gene_biotype = gene_biotype[1],
      .groups = "drop"
    )

  # Calculate TPM
  lengths_kb <- dplyr::pull(new_gene_annot, gene_length_kb)
  names(lengths_kb) <- dplyr::pull(new_gene_annot, !!sym(annotation))
  tpm <- transcripts_per_million(gexp_new, lengths_kb)

  # Final assembly
  new_exp_data <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = gexp_new, tpm = tpm),
    colData = sample_annot,
    rowData = tibble::as_tibble(new_gene_annot)
  )

  return(new_exp_data)
}
