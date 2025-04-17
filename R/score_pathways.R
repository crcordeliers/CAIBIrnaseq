#' Score Pathways Using Gene Expression Data
#'
#' This function scores pathways based on gene expression data using various scoring methods such as GSVA, ssGSEA, plage, and zscore.
#'
#' @param exp_data A `SummarizedExperiment` object containing normalized gene expression data.
#' @param pathways A data frame containing pathways with columns `gene_id`, `gene_symbol`, and `pathway`.
#' @param scoring_method The method used to score pathways. Options are `"gsva"`, `"ssgsea"`, `"plage"`, or `"zscore"`. Default is `"gsva"`.
#' @param verbose Logical. If TRUE, progress messages will be printed. Default is TRUE.
#'
#' @returns A data frame containing the scores for the pathways.
#' @export
#'
score_pathways <- function(exp_data, pathways,
                           scoring_method = "gsva",
                           verbose = TRUE) {
  # Extract normalized expression matrix from the SummarizedExperiment object
  mat <- assays(exp_data)$norm

  # Check the gene annotation in the pathways data
  if (any(rownames(mat) %in% pathways$gene_id)) {
    gene_annot <- "gene_id"
  } else if (any(rownames(mat) %in% pathways$gene_symbol)) {
    gene_annot <- "gene_symbol"
  } else {
    stop("`exp_data` uses unknown gene annotation. Try either using ensembl_gene_id or gene_name/gene_symbol.")
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
