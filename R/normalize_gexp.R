#' Normalize Gene Expression Data
#'
#' This function normalizes gene expression data using either the `rlog` or `vst` method,
#' depending on the number of samples. For fewer than 30 samples, `rlog` normalization is
#' applied; for 30 or more samples, `vst` normalization is used.
#'
#' @param exp_data A `SummarizedExperiment` object containing the gene expression counts
#'        in the `counts` assay.
#'
#' @returns A `SummarizedExperiment` object with the normalized data stored in the `norm` assay.
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
