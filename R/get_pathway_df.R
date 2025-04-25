#' Extract pathway scores and metadata from a SummarizedExperiment
#'
#' This function extracts pathway activity scores for specified pathways and merges them
#' with the sample annotations from a SummarizedExperiment object.
#'
#' @param exp_data A `SummarizedExperiment` object containing sample metadata in `colData`.
#' @param pathway_scores A matrix or data frame with pathway scores (rows = pathways, columns = samples).
#' @param pathways A character vector of pathway names or identifiers to extract.
#'
#' @return A data frame containing the selected pathway scores and sample annotations.
#' @export
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
get_pathway_df <- function(exp_data, pathway_scores, pathways) {
  # --- Input checks ---
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("`exp_data` must be a SummarizedExperiment object.")
  }

  if (!is.character(pathways) || length(pathways) == 0) {
    stop("`pathways` must be a non-empty character vector.")
  }

  if (!is.matrix(pathway_scores) && !is.data.frame(pathway_scores)) {
    stop("`pathway_scores` must be a matrix or data.frame.")
  }

  # --- Sample annotations ---
  sample_annot <- SummarizedExperiment::colData(exp_data) |> as.data.frame()

  # --- Validate pathway presence ---
  missing_pathways <- setdiff(pathways, rownames(pathway_scores))
  if (length(missing_pathways) > 0) {
    warning("The following pathways were not found and will be ignored: ",
            paste(missing_pathways, collapse = ", "))
    pathways <- intersect(pathways, rownames(pathway_scores))
  }

  if (length(pathways) == 0) {
    stop("None of the specified pathways were found in the pathway score matrix.")
  }

  # --- Extract pathway score matrix ---
  if (length(pathways) > 1) {
    path_df <- pathway_scores[pathways, , drop = FALSE] |>
      t() |>
      as.data.frame()
  } else {
    path_vec <- pathway_scores[pathways, ] |> as.numeric()
    path_df <- matrix(
      path_vec, ncol = 1,
      dimnames = list(colnames(pathway_scores), pathways)
    ) |> as.data.frame()
  }

  # --- Merge with sample metadata ---
  path_df <- path_df |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(sample_annot, by = "sample_id")

  return(path_df)
}
