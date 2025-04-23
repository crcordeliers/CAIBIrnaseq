#' Normalize Gene Expression Data
#'
#' This function normalizes gene expression counts in a `SummarizedExperiment` object using either `rlog` or `vst` normalization, depending on the number of samples.
#'
#' @param exp_data A `SummarizedExperiment` object containing the gene expression data with a `counts` assay.
#'
#' @return A `SummarizedExperiment` object with an additional `norm` assay containing the normalized expression data.
#'
#' @details
#' The function applies one of two normalization methods:
#' - `rlog` (regularized log transformation) is used for datasets with fewer than 30 samples.
#' - `vst` (variance-stabilizing transformation) is used for datasets with 30 or more samples.
#'
#' The normalized data is stored in a new assay named `"norm"`.
#'
#' @export
#'
#'@importFrom DESeq2 rlog vst
#'@importFrom SummarizedExperiment assays
normalize_gexp <- function(exp_data) {
  # Extract the counts data from the SummarizedExperiment object
  counts <- SummarizedExperiment::assays(exp_data)[["counts"]]

  # Check if the counts assay is available
  if (is.null(counts)) {
    stop("No 'counts' assay found in the SummarizedExperiment object.")
  }

  # Check the number of samples (columns)
  if (ncol(counts) < 30) {
    message("- Less than 30 samples -> Performing `rlog` normalization...")
    norm <- DESeq2::rlog(counts)  # Apply rlog normalization for small sample sizes
  } else {
    message("- Performing `vst` normalization...")
    norm <- DESeq2::vst(counts)   # Apply vst normalization for large sample sizes
  }

  # Store the normalized data in a new assay in the SummarizedExperiment object
  SummarizedExperiment::assays(exp_data)[["norm"]] <- norm

  # Return the modified SummarizedExperiment object with the normalized data
  return(exp_data)
}
