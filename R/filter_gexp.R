#' Filter Gene Expression Data
#'
#' This function filters genes from a `SummarizedExperiment` object based on minimum expression thresholds.
#'
#' @param exp_data A `SummarizedExperiment` object containing the gene expression data with a `counts` assay.
#' @param min_nsamp An integer specifying the minimum number of samples in which a gene must have expression above `min_counts` to be retained. Default is 1.
#' @param min_counts An integer specifying the minimum count threshold a gene must have in `min_nsamp` samples to be retained. Default is 10.
#'
#' @return A filtered `SummarizedExperiment` object containing only the genes that meet the specified criteria. Two new columns are added to `colData`:
#' - `ncounts`: Total counts per sample.
#' - `nfeats`: Number of features (genes) detected per sample.
#'
#' @details
#' This function removes genes that do not meet the specified thresholds for expression. It adds sample-level metrics (`ncounts` and `nfeats`) to the `colData` of the `SummarizedExperiment` object for downstream analysis.
#'
#' The filtering criteria are:
#' - A gene must have expression greater than or equal to `min_counts` in at least `min_nsamp` samples.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assays colData
filter_gexp <- function(exp_data, min_nsamp = 1, min_counts = 10) {
  # --- VALIDATION CHECKS ---
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("Error: `exp_data` must be a SummarizedExperiment object.")
  }

  counts <- SummarizedExperiment::assays(exp_data)[["counts"]]
  if (is.null(counts)) {
    stop("Error: `counts` assay not found in the SummarizedExperiment object.")
  }

  if (!is.numeric(counts)) {
    stop("Error: `counts` must be a numeric matrix.")
  }

  if (!is.numeric(min_nsamp) || min_nsamp < 1) {
    stop("Error: `min_nsamp` must be a positive integer.")
  }

  if (!is.numeric(min_counts) || min_counts < 0) {
    stop("Error: `min_counts` must be a non-negative integer.")
  }

  # --- FILTERING LOGIC ---
  SummarizedExperiment::colData(exp_data)[["ncounts"]] <- base::colSums(counts)
  SummarizedExperiment::colData(exp_data)[["nfeats"]] <- base::colSums(counts > 0)

  keep_genes <- base::rowSums(counts >= min_counts) >= min_nsamp

  message("- Keeping ", sum(keep_genes), "/", nrow(counts),
          " genes found in at least ", min_nsamp, " sample(s) with at least ",
          min_counts, " counts")

  return(exp_data[keep_genes, ])
}
