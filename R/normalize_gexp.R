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
  counts <- assays(exp_data)[["counts"]]
  if (is.null(counts)) {
    stop("No 'counts' assay found in the SummarizedExperiment object.")
  }
  if (!is.integer(counts)) {
    counts <- matrix(as.integer(round(counts)), nrow = nrow(counts),
                     dimnames = dimnames(counts))
  }
  if (ncol(counts) < 30) {
    message("- Less than 30 samples -> Performing `rlog` normalization...")
    norm <- rlog(counts)
  } else {
    message("- Performing `vst` normalization...")
    norm <- vst(counts)
  }
  assays(exp_data)[["norm"]] <- norm

    return(exp_data)
}
