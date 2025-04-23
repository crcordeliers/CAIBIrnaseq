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
  message("-- Rebasing the gene expression matrix using `", annotation, "` as main annotation")

  # Extract the gene annotations and sample annotations
  gene_annot <- as.data.frame(SummarizedExperiment::rowData(exp_data))
  sample_annot <- SummarizedExperiment::colData(exp_data)

  # Preprocess the gene expression data (assuming 'gexp_preprocess' is defined elsewhere)
  gexp_new <- gexp_preprocess(
    gexp = SummarizedExperiment::assays(exp_data)[["counts"]],
    gene_annotation = gene_annot,
    og_annot = "gene_id",
    keep_annot = annotation,
    keep_stat = "sum"
  )

  gexp_new <- as.matrix(gexp_new)

  # Aggregate annotations based on the specified annotation (e.g., gene_name)
  new_gene_annot <- gene_annot %>%
    dplyr::group_by(!!sym(annotation)) %>%
    dplyr::summarize(gene_id = paste(gene_id, collapse = ", "),
                     gene_length_kb = mean(gene_length_kb),
                     gene_description = paste(gene_description, collapse = ", "),
                     gene_biotype = gene_biotype[1])

  # Calculate TPM (Transcripts Per Million)
  tpm <- transcripts_per_million(gexp_new,
                                 pull(new_gene_annot, gene_length_kb, !!sym(annotation)))

  # Create the new gene annotation dataframe
  gene_annot_df <- tibble::as_tibble(new_gene_annot)

  # Create a new SummarizedExperiment object with the updated data
  new_exp_data <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = gexp_new, tpm = tpm),
    colData = sample_annot,
    rowData = new_gene_annot
  )

  return(new_exp_data)
}
