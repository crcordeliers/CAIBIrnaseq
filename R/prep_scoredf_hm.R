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
  available_paths <- rownames(pathway_scores)

  if (all(is.numeric(pathways))) {
    valid_indices <- intersect(seq_along(available_paths), pathways)
    feats <- available_paths[valid_indices]
  } else if (all(pathways %in% available_paths)) {
    feats <- pathways
  } else {
    stop("Not all specified pathways are available in the provided pathway scores.")
  }

  path_table <- pathway_scores[feats, , drop = FALSE] |>
    t() |>
    as.data.frame()

  # Remove duplicate 'sample_id' if already present
  path_table <- path_table[, setdiff(names(path_table), "sample_id")]

  path_table <- tibble::rownames_to_column(path_table, "sample_id")


  return(list(
    table = path_table,
    colv = "sample_id",
    rowv = feats
  ))
}
