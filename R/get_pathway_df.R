#' Extract pathway scores and metadata from a SummarizedExperiment
#'
#' Extract pathway activity scores and sample annotations
#'
#' This function extracts pathway scores for specified pathways and merges them
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
#'
get_pathway_df <- function(exp_data, pathway_scores, pathways) {
  sample_annot <- SummarizedExperiment::colData(exp_data) |> as.data.frame()

  if (length(pathways) > 1) {
    path_df <- pathway_scores[pathways, ] |>
      t() |>
      as.data.frame()
  } else {
    path_vec <- pathway_scores[pathways, ] |> as.numeric()
    path_df <- matrix(
      path_vec, ncol = 1,
      dimnames = list(colnames(pathway_scores), pathways)
    ) |> as.data.frame()
  }

  path_df <- path_df |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(sample_annot, by = "sample_id")

  return(path_df)
}
