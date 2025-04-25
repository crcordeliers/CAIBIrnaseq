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
  # Basic checks
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("`exp_data` must be a 'SummarizedExperiment' object.")
  }

  if (!is.matrix(pathway_scores) && !is.data.frame(pathway_scores)) {
    stop("`pathway_scores` must be a matrix or data frame.")
  }

  if (is.null(rownames(pathway_scores))) {
    stop("`pathway_scores` must have row names representing pathway names.")
  }

  if (!is.numeric(pathways) && !is.character(pathways)) {
    stop("`pathways` must be either numeric indices or character pathway names.")
  }

  available_paths <- rownames(pathway_scores)

  # Validate and select pathways
  if (is.numeric(pathways)) {
    valid_indices <- intersect(seq_along(available_paths), pathways)
    if (length(valid_indices) == 0) {
      stop("None of the specified indices are valid for `pathway_scores`.")
    }
    feats <- available_paths[valid_indices]
  } else if (all(pathways %in% available_paths)) {
    feats <- pathways
  } else {
    stop("Some of the specified pathway names are not found in `pathway_scores`.")
  }

  # Extract and clean sample annotations
  samp_annot <- SummarizedExperiment::colData(exp_data) |> as.data.frame()

  if ("sample_id" %in% colnames(samp_annot)) {
    warning("'sample_id' already exists in colData and will be overwritten.")
    samp_annot <- samp_annot[, colnames(samp_annot) != "sample_id"]
  }

  samp_annot <- tibble::rownames_to_column(samp_annot, "sample_id")

  # Prepare tidy pathway score table
  path_table <- pathway_scores[feats, , drop = FALSE] |>
    t() |>
    as.data.frame()

  if ("sample_id" %in% colnames(path_table)) {
    warning("'sample_id' already exists in `pathway_scores` and will be overwritten.")
    path_table <- path_table[, colnames(path_table) != "sample_id"]
  }

  path_table <- tibble::rownames_to_column(path_table, "sample_id")

  # Join with sample annotations
  path_table <- dplyr::left_join(path_table, samp_annot, by = "sample_id")

  return(list(
    table = path_table,
    colv = "sample_id",
    rowv = feats
  ))
}
