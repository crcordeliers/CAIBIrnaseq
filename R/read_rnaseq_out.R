#' Read RNA-seq output from nf-core/rnaseq
#'
#' This function loads and formats a `SummarizedExperiment` object
#' from the nf-core/rnaseq output directory (e.g., generated using STAR + Salmon).
#' It sets up TPM, log-transformed TPM, cleans the sample annotations,
#' and adds basic QC metadata (number of counts and detected features).
#'
#' @param DATA_PATH A character string pointing to the path containing nf-core/rnaseq results.
#'
#' @return A `SummarizedExperiment` object with assays and annotations formatted.
#' @export
#'
#' @importFrom readr read_rds
#' @importFrom SummarizedExperiment assays colData rowData
#'
read_rnaseq_out <- function(DATA_PATH) {
  # Load the rds file from nf-core/rnaseq
  exp_data <- readr::read_rds(file.path(DATA_PATH, "star_salmon/salmon.merged.gene_counts_length_scaled.rds"))

  # Rename second assay to 'tpm'
  names(SummarizedExperiment::assays(exp_data))[2] <- "tpm"

  # Add log-transformed TPM
  SummarizedExperiment::assays(exp_data)$log_tpm <- log2(SummarizedExperiment::assays(exp_data)$tpm + 1)

  # Clean up colData
  if ("files" %in% colnames(colData(exp_data))) {
    SummarizedExperiment::colData(exp_data)$files <- NULL
  }
  colnames(SummarizedExperiment::colData(exp_data)) <- "sample_id"

  # Round counts to integers
  SummarizedExperiment::assays(exp_data)[["counts"]] <- as.matrix(
    round(SummarizedExperiment::assays(exp_data)[["counts"]])
  )

  # Add QC metadata
  counts <- SummarizedExperiment::assays(exp_data)[["counts"]]
  SummarizedExperiment::rowData(exp_data)$ncounts <- base::rowSums(counts)
  SummarizedExperiment::colData(exp_data)$nfeats <- apply(counts, 2, function(x) sum(x > 0))
  SummarizedExperiment::colData(exp_data)$ncounts <- base::colSums(counts)

  return(exp_data)
}
