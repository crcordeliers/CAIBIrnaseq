#' Differential Expression Analysis
#'
#' Performs differential expression analysis using DESeq2 on a given count matrix and sample information.
#'
#' @param countData A matrix or data frame of raw count data. Rows represent genes, and columns represent samples.
#' @param sampleInfo A data frame containing sample metadata. Must include a `condition` column specifying the experimental conditions.
#' @param method A string specifying the method for differential expression analysis. Currently supports only `"DESeq2"`. Default is `"DESeq2"`.
#' @param cutoff An integer specifying the minimum number of counts required across all samples for a gene to be included in the analysis. Default is `10`.
#' @param design A
#' @param coefname A
#'
#' @details
#' This function performs differential expression analysis using the DESeq2 package. It filters genes with low counts, estimates size factors for normalization, and performs the DESeq2 analysis pipeline. Log fold-change shrinkage is applied using the `lfcShrink` function.
#'
#' @return A data frame containing the results of the differential expression analysis, including adjusted p-values, log fold changes, and other statistics.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink
#' @importFrom BiocGenerics estimateSizeFactors counts
#' @importFrom SummarizedExperiment colData
#'
#' @export
diffExpAnalysis <- function(countData, sampleInfo, method = "DESeq2", cutoff = 10, design, coefname) {
  # DESeq2 Method
  if (tolower(method) == "deseq2") {
    # Create DESeqDataSet object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                          colData = sampleInfo,
                                          design = ~ Genotype)

    keep <- rowSums(BiocGenerics::counts(dds)) >= cutoff
    dds <- dds[keep,]
    dds <- BiocGenerics::estimateSizeFactors(dds)

    # Perform DESeq2 analysis
    dds <- DESeq2::DESeq(dds)

    print(DESeq2::resultsNames(dds))

    res <- DESeq2::results(dds)
    res_shrink <- DESeq2::lfcShrink(dds, coef = coefname, type = "normal")

    return(as.data.frame(res_shrink))
  }
}
