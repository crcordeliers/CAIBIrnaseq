#' Get Pathway Genes
#'
#' Extracts genes associated with specific pathways from an annotation collection.
#'
#' @param annotation_collection A data frame containing pathway annotations. It should have at least two columns: `pathway` and a gene annotation column (e.g., `gene_symbol` or `gene_id`).
#' @param pathway_names A character vector specifying the names of pathways for which to retrieve associated genes.
#' @param annotation A character string specifying the column in `annotation_collection` to use for gene annotations. Default is `"gene_symbol"`.
#'
#' @details
#' This function filters the provided annotation collection for the specified pathway names and extracts the corresponding genes based on the selected annotation column.
#'
#' @return A character vector of genes associated with the specified pathways.
#'
#' @importFrom dplyr filter pull
#' @importFrom rlang sym
#' @export
get_pathway_genes <- function(annotation_collection,
                              pathway_names,
                              annotation = "gene_symbol") {
  # --- Input checks ---
  if (!is.data.frame(annotation_collection)) {
    stop("`annotation_collection` must be a data frame.")
  }

  if (!"pathway" %in% colnames(annotation_collection)) {
    stop("`annotation_collection` must contain a `pathway` column.")
  }

  if (!annotation %in% colnames(annotation_collection)) {
    stop("The column '", annotation, "' was not found in `annotation_collection`.")
  }

  if (!is.character(pathway_names) || length(pathway_names) == 0) {
    stop("`pathway_names` must be a non-empty character vector.")
  }

  # --- Check for missing pathways ---
  missing_paths <- setdiff(pathway_names, unique(annotation_collection$pathway))
  if (length(missing_paths) > 0) {
    warning("Some pathways were not found in the annotation collection and will be ignored: ",
            paste(missing_paths, collapse = ", "))
  }

  # --- Extract genes ---
  genes <- annotation_collection |>
    dplyr::filter(.data$pathway %in% pathway_names) |>
    dplyr::pull(!!rlang::sym(annotation))

  return(unique(genes))
}
