#' Differential Expression Analysis
#'
#' Performs differential expression analysis using DESeq2 on a given count matrix and sample information.
#'
#' @param countData A matrix or data frame of raw count data. Rows represent genes, and columns represent samples.
#' @param sampleInfo A data frame containing sample metadata. Must include a `condition` column specifying the experimental conditions.
#' @param method A string specifying the method for differential expression analysis. Currently supports only `"DESeq2"`. Default is `"DESeq2"`.
#' @param cutoff An integer specifying the minimum number of counts required across all samples for a gene to be included in the analysis. Default is `10`.
#' @param annotation ee

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
diffExpAnalysis <- function(countData, sampleInfo, method = "DESeq2", cutoff = 10, annotation) {
  # DESeq2 Method
  if (tolower(method) == "deseq2") {
    # Create DESeqDataSet object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                          colData = sampleInfo,
                                  design = as.formula(paste('~',annotation)))

    keep <- rowSums(BiocGenerics::counts(dds)) >= cutoff
    dds <- dds[keep,]

    dds <- BiocGenerics::estimateSizeFactors(dds)

    # Perform DESeq2 analysis
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    res_shrink <- DESeq2::lfcShrink(dds, coef = tail(DESeq2::resultsNames(dds), 1))

    return(as.data.frame(res_shrink))

  }
  # Limma-Voom Method
  else if (tolower(method) == "limma") {
    # Create DGEList object for Limma
    dge <- DGEList(counts = countData)
    dge <- calcNormFactors(dge)
    drop <- which(apply(cpm(dge), 1, max) < cutoff)
    dge <- dge[-drop,]

    v <- voom(dge, design = model.matrix(~ condition, data = sampleInfo), plot = TRUE)

    fit <- lmFit(v, model.matrix(~ condition, data = sampleInfo))

    fit <- eBayes(fit)

    top_table <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

    return(top_table)

  }
  else {
    stop("Invalid method. Choose either 'DESeq2' or 'Limma'.")
  }
}
