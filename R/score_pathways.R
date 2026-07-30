#' Score Pathways
#'
#' Computes pathway scores for a given expression dataset using specified scoring methods.
#'
#' @param exp_data A `SummarizedExperiment` object containing normalized expression data in the `assays(exp_data)$norm` slot.
#' @param pathways Either a data frame with pathway definitions, containing at least two columns:
#' `pathway` (pathway name) and either `gene_id` (Ensembl IDs) or `gene_symbol` (gene symbols); or
#' a named list of character vectors (one gene-ID vector per pathway/signature), e.g.
#' `list(MySignature = c("GENE1", "GENE2"))`. The list form lets custom gene signatures be scored
#' the same way as an MSigDB collection - the gene IDs must match the rownames of
#' `assays(exp_data)$norm` directly (no gene_id/gene_symbol auto-detection is performed).
#' @param scoring_method A character string specifying the scoring method to use. Options are `"gsva"`, `"ssgsea"`, `"plage"`, or `"zscore"`. Default is `"gsva"`.
#' @param min_genes An integer specifying the minimum number of genes required for a pathway to be considered. Default is `100`. Lower only if using targetted panel.
#' @param verbose Logical; if `TRUE`, prints progress messages during computation. Default is `TRUE`.
#'
#' @details
#' When `pathways` is a data frame, the function identifies the gene annotation used in the
#' expression matrix (`gene_id` or `gene_symbol`) by matching row names of `assays(exp_data)$norm`
#' to the `pathways` data frame. When `pathways` is a named list, its gene IDs are used as-is
#' (matched directly against `rownames(assays(exp_data)$norm)`), which is the natural format for
#' scoring custom signatures rather than an MSigDB collection.
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
#' @importFrom matrixStats rowVars
#'
#' @export
score_pathways <- function(exp_data, pathways,
                           scoring_method = "gsva",
                           min_genes = 100,
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

  if (is.data.frame(pathways)) {
    # Check the gene annotation in the pathways data
    votes_gene_id <- sum(rownames(mat) %in% pathways$gene_id)
    votes_gene_symbol <- sum(rownames(mat) %in% pathways$gene_symbol)

    if (votes_gene_id >= votes_gene_symbol & votes_gene_id > min_genes) {
      gene_annot <- "gene_id"
    } else if (votes_gene_symbol > votes_gene_id & votes_gene_symbol > min_genes) {
      gene_annot <- "gene_symbol"
    } else {
      stop("`exp_data` uses unknown gene annotation.
           Try either using ensembl_gene_id or gene_name/gene_symbol.")
    }

    gene_sets <- .as_gene_sets(pathways, id_col = gene_annot)
  } else {
    # Named list of custom gene signatures: use the gene IDs as-is.
    gene_sets <- .as_gene_sets(pathways)
    matched <- sum(rownames(mat) %in% unique(unlist(gene_sets, use.names = FALSE)))
    if (matched < min_genes) {
      stop("`exp_data` uses unknown gene annotation: fewer than `min_genes` (", min_genes, ") of ",
           "`rownames(assays(exp_data)$norm)` match the gene IDs in `pathways`. Make sure they use ",
           "the same gene ID type (e.g. both Ensembl IDs, or both gene symbols).")
    }
  }

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

  # Sort pathways by most variable scores
  path_vars <- rowVars(as.matrix(path_scores_df))
  names(path_vars) <- rownames(path_scores_df)
  max_diff_paths <- path_vars |>
    sort(decreasing = TRUE) |>
    names()

  # Filter the data frame to only include the most variable pathways
  path_scores_df <- path_scores_df[max_diff_paths, ]

  return(path_scores_df)
}
