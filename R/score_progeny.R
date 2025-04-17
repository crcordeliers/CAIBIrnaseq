#' Score PROGENy Pathways for Gene Expression Data
#'
#' This function calculates pathway scores for gene expression data using the PROGENy algorithm.
#' PROGENy is an algorithm that predicts pathway activity based on gene expression data.
#'
#' @param exp_data A `SummarizedExperiment` object containing normalized gene expression data.
#' @param species The species for which to score the pathways. Defaults to "Homo sapiens" ("Human").
#'
#' @returns A data frame containing the PROGENy pathway scores for each sample.
#' @export
#'
#' @importFrom progeny progeny
#'
score_progeny <- function(exp_data, species = "Homo sapiens") {
  # Extract the normalized gene expression matrix
  gexp <- assays(exp_data)[["norm"]]

  # Determine the organism type based on species input
  organism <- if_else(species == "Homo sapiens", "Human", "Mouse")

  # Calculate PROGENy pathway scores
  progeny_scores <- t(progeny(gexp, organism = organism))  # transposing for correct format

  # Convert the resulting scores into a data frame
  progeny_scores_df <- data.frame(progeny_scores)

  # Return the pathway scores as a data frame
  return(progeny_scores_df)
}
