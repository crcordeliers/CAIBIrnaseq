#' Differential Expression Analysis
#'
#' Performs differential expression analysis using DESeq2 on a given count matrix and sample information.
#'
#' @param countData A matrix or data frame of raw count data. Rows represent genes, and columns represent samples.
#' @param sampleInfo A data frame containing sample metadata. Must include a `condition` column specifying the experimental conditions.
#' @param method A string specifying the method for differential expression analysis. Currently supports only `"DESeq2"`. Default is `"DESeq2"`.
#' @param cutoff An integer specifying the minimum number of counts required across all samples for a gene to be included in the analysis. Default is `10`.
#' @param design A formula specifying the experimental design (e.g., ~ condition).
#' @param coefname A string specifying the coefficient name to shrink. Must match one of the results names in the model.
#'
#' @details
#' This function performs differential expression analysis using the DESeq2 package. It filters genes with low counts, estimates size factors for normalization, and performs the DESeq2 analysis pipeline. Log fold-change shrinkage is applied using the `lfcShrink` function.
#'
#' @return A data frame containing the results of the differential expression analysis, including adjusted p-values, log fold changes, and other statistics.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink resultsNames
#' @importFrom BiocGenerics estimateSizeFactors counts
#' @importFrom SummarizedExperiment colData
#'
#' @export
diffExpAnalysis <- function(countData, sampleInfo, method = "DESeq2", cutoff = 10, design, coefname) {

  # ----- VALIDATION CHECKS -----

  # Check method
  if (tolower(method) != "deseq2") {
    stop("Error: Only 'DESeq2' is currently supported.")
  }

  # Check countData
  if (!is.matrix(countData) && !is.data.frame(countData)) {
    stop("Error: `countData` must be a matrix or data frame.")
  }

  # Check sampleInfo
  if (!is.data.frame(sampleInfo)) {
    stop("Error: `sampleInfo` must be a data frame.")
  }

  # Check that rownames(sampleInfo) match colnames(countData)
  if (!all(colnames(countData) %in% rownames(sampleInfo))) {
    stop("Error: Column names of `countData` must match row names of `sampleInfo`.")
  }

  # Check cutoff
  if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff < 0) {
    stop("Error: `cutoff` must be a non-negative numeric value.")
  }

  # Check coefname
  if (missing(coefname) || !is.character(coefname) || length(coefname) != 1) {
    stop("Error: `coefname` must be a single character string.")
  }

  # ----- ANALYSIS -----
  # Create DESeq2 dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData = sampleInfo,
    design = design
  )

  # Filter low-expressed genes
  keep <- rowSums(BiocGenerics::counts(dds)) >= cutoff
  dds <- dds[keep, ]

  # Normalize
  dds <- BiocGenerics::estimateSizeFactors(dds)

  # Run DESeq
  dds <- DESeq2::DESeq(dds)

  # Print available coefficient names
  available_coefs <- DESeq2::resultsNames(dds)
  print(available_coefs)

  # Check if coefname exists
  if (!(coefname %in% available_coefs)) {
    stop("Error: The specified `coefname` does not match any available coefficients.\nAvailable coefficients are:\n",
         paste(available_coefs, collapse = ", "))
  }

  # Log fold-change shrinkage
  res_shrink <- DESeq2::lfcShrink(dds, coef = coefname, type = "normal")

  return(as.data.frame(res_shrink))
}
