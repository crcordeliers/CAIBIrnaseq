#' Score PROGENy Pathways
#'
#' Calculate PROGENy pathway activity scores based on normalized gene expression data.
#'
#' @param exp_data A `SummarizedExperiment` object containing normalized gene expression data.
#' @param species A character string specifying the species, either `"Homo sapiens"` (default) or `"Mus musculus"`.
#'
#' @return A data frame containing PROGENy pathway activity scores for each sample.
#'
#' @details
#' The function uses the normalized gene expression data from the specified assay of the `SummarizedExperiment` object to compute PROGENy scores.
#' PROGENy scores represent pathway activity and are computed based on predefined pathway models for the specified species.
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom progeny progeny
#' @importFrom dplyr if_else
#'
#' @export
#'
score_progeny <- function(exp_data, species = "Homo sapiens") {
  # Extract the normalized gene expression matrix
  gexp <- SummarizedExperiment::assays(exp_data)[["norm"]]

  # Determine the organism type based on species input
  organism <- dplyr::if_else(species == "Homo sapiens", "Human", "Mouse")

  # Calculate PROGENy pathway scores
  progeny_scores <- t(progeny::progeny(gexp, organism = organism))  # transposing for correct format

  # Convert the resulting scores into a data frame
  progeny_scores_df <- data.frame(progeny_scores)

  # Return the pathway scores as a data frame
  return(progeny_scores_df)
}
