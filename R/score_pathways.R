#' Score Pathways
#'
#' Computes pathway scores for a given expression dataset using specified scoring methods.
#'
#' @param exp_data A `SummarizedExperiment` object containing normalized expression data in the `assays(exp_data)$norm` slot.
#' @param pathways A data frame with pathway definitions, containing at least two columns: `pathway` (pathway name) and either `gene_id` (Ensembl IDs) or `gene_symbol` (gene symbols).
#' @param scoring_method A character string specifying the scoring method to use. Options are `"gsva"`, `"ssgsea"`, `"plage"`, or `"zscore"`. Default is `"gsva"`.
#' @param verbose Logical; if `TRUE`, prints progress messages during computation. Default is `TRUE`.
#'
#' @details
#' The function identifies the gene annotation used in the expression matrix (`gene_id` or `gene_symbol`) by matching row names of `assays(exp_data)$norm` to the `pathways` data frame.
#' It then splits the pathways into gene sets and scores them using the specified method from the `GSVA` package.
#'
#' The available scoring methods are:
#' \describe{
#'   \item{`gsva`}{Gene Set Variation Analysis.}
#'   \item{`ssgsea`}{Single-sample Gene Set Enrichment Analysis.}
#'   \item{`plage`}{Pathway Level Analysis of Gene Expression.}
#'   \item{`zscore`}{Z-score normalization.}
#' }
#'
#' Pathway scores are sorted by the sum of absolute scores across samples, prioritizing pathways with the highest variation.
#'
#' @return A data frame of pathway scores, with pathways as row names and samples as columns. Pathways are sorted by their total score variation.
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom GSVA gsvaParam ssgseaParam plageParam zscoreParam gsva
#'
#' @export
score_pathways <- function(exp_data, pathways,
                           scoring_method = "gsva",
                           verbose = TRUE) {


  # Check if the input data is of the correct class
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("`exp_data` must be a SummarizedExperiment object.")
  }

  # Check if the necessary assay exists
  if (!"norm" %in% names(assays(exp_data))) {
    stop("`exp_data` does not contain the normalized expression data in `assays(exp_data)$norm`.")
  }

  # Extract normalized expression matrix from the SummarizedExperiment object
  mat <- assays(exp_data)$norm

  # Check the gene annotation in the pathways data
  if (any(rownames(mat) %in% pathways$gene_id)) {
    gene_annot <- "gene_id"
  } else if (any(rownames(mat) %in% pathways$gene_symbol)) {
    gene_annot <- "gene_symbol"
  } else {
    stop("`exp_data` uses unknown gene annotation.
         Try either using ensembl_gene_id or gene_name/gene_symbol.")
  }

  # Generate a list of gene sets based on the pathway annotations
  gene_sets <- split(pathways[[gene_annot]], pathways$pathway)

  # Score pathways using the specified method
  if (scoring_method == "gsva") {
    param <- GSVA::gsvaParam(mat, geneSets = gene_sets, minSize = 2)
  } else if (scoring_method == "ssgsea") {
    param <- GSVA::ssgseaParam(mat, geneSets = gene_sets, minSize = 2)
  } else if (scoring_method == "plage") {
    param <- GSVA::plageParam(mat, geneSets = gene_sets, minSize = 2)
  } else if (scoring_method == "zscore") {
    param <- GSVA::zscoreParam(mat, geneSets = gene_sets, minSize = 2)
  } else {
    stop("Invalid scoring method. Choose from 'gsva', 'ssgsea', 'plage', or 'zscore'.")
  }

  # Perform the GSVA analysis to get pathway scores
  scores <- GSVA::gsva(param, verbose = verbose)
  path_scores_df <- data.frame(scores)

  # Sort pathways by absolute scores and select the most variable pathways
  max_diff_paths <- path_scores_df |>
    abs() |>
    rowSums() |>
    sort(decreasing = TRUE) |>
    names()

  # Filter the data frame to only include the most variable pathways
  path_scores_df <- path_scores_df[max_diff_paths, ]

  return(path_scores_df)
}
