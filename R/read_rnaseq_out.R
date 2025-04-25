#' Read RNA-seq Output from nf-core/rnaseq
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
  if (!is.character(DATA_PATH) || length(DATA_PATH) != 1) {
    stop("`DATA_PATH` must be a single character string.")
  }

  rds_path <- file.path(DATA_PATH, "star_salmon/salmon.merged.gene_counts_length_scaled.rds")
  if (!file.exists(rds_path)) {
    stop("The RDS file was not found at the expected path: ", rds_path)
  }

  # Load SummarizedExperiment
  exp_data <- readr::read_rds(rds_path)

  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("The loaded object is not a `SummarizedExperiment`.")
  }

  # Add log-transformed TPM
  tpm <- SummarizedExperiment::assays(exp_data)[["tpm"]]
  SummarizedExperiment::assays(exp_data)[["log_tpm"]] <- log2(tpm + 1)

  # Clean up colData: remove 'files' column if present and rename all columns to 'sample_id' if appropriate
  cd <- SummarizedExperiment::colData(exp_data)

  if ("files" %in% colnames(cd)) {
    cd$files <- NULL
  }

  # Rename only if there's a single column — avoid overwriting existing metadata
  if (ncol(cd) == 1) {
    colnames(cd) <- "sample_id"
  }

  SummarizedExperiment::colData(exp_data) <- cd

  # Round counts to integers
  counts <- SummarizedExperiment::assays(exp_data)[["counts"]]
  if (!is.null(counts)) {
    SummarizedExperiment::assays(exp_data)[["counts"]] <- as.matrix(round(counts))
  } else {
    stop("`counts` assay is missing from the object.")
  }

  # Add QC metadata
  SummarizedExperiment::rowData(exp_data)$ncounts <- rowSums(counts)
  SummarizedExperiment::colData(exp_data)$nfeats <- colSums(counts > 0)
  SummarizedExperiment::colData(exp_data)$ncounts <- colSums(counts)

  return(exp_data)
}
