#' Prepare Pathway Scores Data Frame for Heatmap
#'
#' This function prepares a data frame of pathway scores, transposed and ready for
#' heatmap visualization. It filters and selects a subset of pathways, either by
#' name or index.
#'
#' @param pathway_scores A matrix or data frame with pathway scores (rows = pathways, columns = samples).
#' @param pathways Either a numeric vector of row indices or a character vector of pathway names to include. Default is 1:20.
#'
#' @returns A list with the following elements:
#' \describe{
#'   \item{table}{A data frame in wide format (samples in rows, pathways in columns) with `sample_id` as row identifier.}
#'   \item{colv}{The column name to use as sample ID (`"sample_id"`).}
#'   \item{rowv}{A character vector of selected pathway names.}
#' }
#' @export
#'
#' @importFrom tibble rownames_to_column
#'
prep_scoredf_hm <- function(pathway_scores, pathways = 1:20) {
  # Basic checks
  if (!is.matrix(pathway_scores) && !is.data.frame(pathway_scores)) {
    stop("`pathway_scores` must be a matrix or a data frame.")
  }

  if (is.null(rownames(pathway_scores))) {
    stop("`pathway_scores` must have row names representing pathway names.")
  }

  if (!is.numeric(pathways) && !is.character(pathways)) {
    stop("`pathways` must be either numeric (row indices) or character (pathway names).")
  }

  available_paths <- rownames(pathway_scores)

  if (is.numeric(pathways)) {
    valid_indices <- intersect(seq_along(available_paths), pathways)
    if (length(valid_indices) == 0) {
      stop("None of the provided indices are valid for `pathway_scores`.")
    }
    feats <- available_paths[valid_indices]
  } else if (all(pathways %in% available_paths)) {
    feats <- pathways
  } else {
    stop("Not all specified pathways are available in the provided `pathway_scores`.")
  }

  # Subset, transpose, and format
  path_table <- pathway_scores[feats, , drop = FALSE] |>
    t() |>
    as.data.frame()

  # Avoid duplicate sample_id column
  path_table <- path_table[, setdiff(names(path_table), "sample_id")]

  path_table <- tibble::rownames_to_column(path_table, "sample_id")

  return(list(
    table = path_table,
    colv = "sample_id",
    rowv = feats
  ))
}
