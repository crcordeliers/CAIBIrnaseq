#' Score PROGENy Pathways
#'
#' Calculate PROGENy pathway activity scores based on normalized gene expression data.
#'
#' @param exp_data A `SummarizedExperiment` object containing normalized gene expression data in the `assays(exp_data)$norm` slot.
#' @param species A character string specifying the species, either `"Homo sapiens"` (default) or `"Mus musculus"`.
#' @param min_genes An integer specifying the minimum number of gene symbols in `exp_data` that must
#' match PROGENy's model genes before scoring proceeds. Default is `100`.
#'
#' @return A data frame containing PROGENy pathway activity scores for each sample.
#'
#' @details
#' The function uses the normalized gene expression data from the specified assay of the `SummarizedExperiment` object to compute PROGENy scores.
#' PROGENy scores represent pathway activity and are computed based on predefined pathway models for the specified species.
#'
#' PROGENy's models are keyed on gene symbols (HGNC for human, MGI for mouse), so `exp_data` must have
#' gene symbols as rownames of the `norm` assay, not Ensembl (or other) gene IDs. Use `rebase_gexp()`
#' beforehand to convert if needed.
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom progeny progeny getModel
#' @importFrom dplyr if_else
#'
#' @export
score_progeny <- function(exp_data, species = "Homo sapiens", min_genes = 100) {
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

  # Validate that `exp_data` is annotated with gene symbols matching PROGENy's model.
  # PROGENy's model matrices are keyed on gene symbols, so Ensembl (or other) IDs would
  # silently fail to match, producing near-empty/NA scores with no error.
  model_genes <- rownames(progeny::getModel(organism = organism))
  n_matched <- sum(rownames(gexp) %in% model_genes)

  if (n_matched < min_genes) {
    looks_ensembl <- mean(grepl("^ENS", rownames(gexp))) > 0.5
    hint <- if (looks_ensembl) {
      " `exp_data` looks like it uses Ensembl gene IDs; convert to gene symbols first with `rebase_gexp()`."
    } else {
      " Make sure `exp_data` uses gene symbols (e.g. via `rebase_gexp()`), not another gene ID type."
    }
    stop(
      "`exp_data` does not appear to use gene symbols matching the PROGENy `", organism, "` model ",
      "(only ", n_matched, " overlapping genes found; need at least ", min_genes, ").", hint
    )
  }

  # Calculate PROGENy pathway scores
  progeny_scores <- t(progeny::progeny(gexp, organism = organism))  # Transpose for correct format

  progeny_scores_df <- data.frame(progeny_scores)

  # Return the pathway scores as a data frame
  return(progeny_scores_df)
}
