#' Filter gene expression data based on counts and sample presence
#'
#' This function filters genes from a `SummarizedExperiment` object based on the
#' number of counts in at least a certain number of samples. Genes that do not
#' meet the specified threshold for both the minimum number of samples and the
#' minimum number of counts per sample are removed.
#'
#' @param exp_data A `SummarizedExperiment` object containing the gene expression data.
#' @param min_nsamp An integer. The minimum number of samples in which a gene must have
#'                  at least `min_counts` to be retained. Default is 1.
#' @param min_counts An integer. The minimum number of counts for a gene to be retained
#'                   in each sample. Default is 10.
#'
#' @returns A `SummarizedExperiment` object containing only the genes that meet the criteria.
#' @export
#'
#' @importFrom SummarizedExperiment colData
#'
#'
filter_gexp <- function(exp_data, min_nsamp = 1, min_counts = 10) {
  # Extract counts matrix from the SummarizedExperiment object
  counts <- assays(exp_data)[["counts"]]

  # Calculate the total counts per sample and the number of non-zero features
  SummarizedExperiment::colData(exp_data)[["ncounts"]] <- base::colSums(counts)
  SummarizedExperiment::colData(exp_data)[["nfeats"]] <- base::colSums(counts > 0)

  # Identify genes to keep based on the filter criteria
  keep_genes <- base::rowSums(counts >= min_counts) >= min_nsamp

  # Print a message showing how many genes are kept
  message("- Keeping ", sum(keep_genes), "/", nrow(counts),
          " genes found in at least ", min_nsamp, " sample(s) with at least ",
          min_counts, " counts")

  # Return the filtered SummarizedExperiment object with the selected genes
  return(exp_data[keep_genes, ])
}
