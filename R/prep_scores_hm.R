#' Prepare Pathway Scores for Heatmap Visualization
#'
#' This function prepares a tidy table of pathway scores joined with sample annotations
#' from a `SummarizedExperiment` object. It returns a list suitable for heatmap plotting.
#'
#' @param exp_data A `SummarizedExperiment` object containing sample annotations in `colData`.
#' @param pathway_scores A matrix or data frame with pathway scores (rows = pathways, columns = samples).
#' @param pathways Either a numeric vector indicating indices of pathways to include, or a character vector of pathway names. Default is 1:20.
#'
#' @returns A list containing:
#' \describe{
#'   \item{table}{A data frame in wide format with pathway scores and sample metadata.}
#'   \item{colv}{Column to use as identifier for samples (usually `"sample_id"`).}
#'   \item{rowv}{A character vector of selected pathways.}
#' }
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment colData
#'
prep_scores_hm <- function(exp_data, pathway_scores, pathways = 1:20) {
  samp_annot <- SummarizedExperiment::colData(exp_data) |>
    as.data.frame()

  # Check if 'sample_id' already exists in colData
  if ("sample_id" %in% colnames(samp_annot)) {
    warning("'sample_id' already exists in colData. It will be overwritten by rownames.")
    samp_annot <- samp_annot[, colnames(samp_annot) != "sample_id"]
  }

  # Convert rownames of sample annotations to 'sample_id'
  samp_annot <- tibble::rownames_to_column(samp_annot, "sample_id")

  available_paths <- rownames(pathway_scores)

  # Validate and subset pathways
  if (all(is.numeric(pathways))) {
    valid_indices <- intersect(1:length(available_paths), pathways)
    feats <- available_paths[valid_indices]
  } else if (all(pathways %in% available_paths)) {
    feats <- pathways
  } else {
    stop("Not all pathways are available in the provided pathway scores.")
  }

  # Create tidy table from pathway scores
  path_table <- pathway_scores[feats, , drop = FALSE] |>
    t() |>
    as.data.frame()

  # Check if 'sample_id' already exists in pathway scores table
  if ("sample_id" %in% colnames(path_table)) {
    warning("'sample_id' already exists in pathway_scores. It will be overwritten by rownames.")
    path_table <- path_table[, colnames(path_table) != "sample_id"]
  }

  # Add 'sample_id' to pathway scores table
  path_table <- tibble::rownames_to_column(path_table, "sample_id")

  # Join pathway scores with sample annotations
  path_table <- dplyr::left_join(path_table, samp_annot, by = "sample_id")

  return(list(table = path_table, colv = "sample_id", rowv = feats))
}
