#' Score PROGENy Pathways
#'
#' Calculate PROGENy pathway activity scores based on normalized gene expression data.
#'
#' @param exp_data A `SummarizedExperiment` object containing normalized gene expression data in the `assays(exp_data)$norm` slot.
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
score_progeny <- function(exp_data, species = "Homo sapiens") {
  # Check if the input is a SummarizedExperiment object
  if (!inherits(exp_data, "SummarizedExperiment")) {
    stop("`exp_data` must be a SummarizedExperiment object.")
  }

  # Check if the 'norm' assay exists in the exp_data
  if (!"norm" %in% names(SummarizedExperiment::assays(exp_data))) {
    stop("`exp_data` does not contain the 'norm' assay with normalized expression data.")
  }

  # Extract the normalized gene expression matrix
  gexp <- SummarizedExperiment::assays(exp_data)[["norm"]]

  if (species == "Homo sapiens" || species == " homo sapiens" || species == "human" || species == "Human") {
    organism <- "Human"
  } else if (species == "Mus musculus" || species == "mus musculus" || species == "Mouse" || species == "mouse") {
    organism <- "Mouse"
  } else {
    stop("Invalid species. Choose either 'Homo sapiens' or 'Mus musculus'.")
  }

  # Calculate PROGENy pathway scores
  progeny_scores <- t(progeny::progeny(gexp, organism = organism))  # Transpose for correct format

  # Convert the resulting scores into a data frame
  progeny_scores_df <- data.frame(progeny_scores)

  # Return the pathway scores as a data frame
  return(progeny_scores_df)
}
