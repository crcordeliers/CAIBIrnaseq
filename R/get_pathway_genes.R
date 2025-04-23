#' Get Pathway Genes
#'
#' Extracts genes associated with specific pathways from an annotation collection.
#'
#' @param annotation_collection A data frame containing pathway annotations. It should have at least two columns: `pathway` (pathway names) and a column for gene annotations (e.g., `gene_symbol` or `gene_id`).
#' @param pathway_names A character vector specifying the names of pathways for which to retrieve associated genes.
#' @param annotation A character string specifying the column in `annotation_collection` to use for gene annotations. Default is `"gene_symbol"`.
#'
#' @details
#' This function filters the provided annotation collection for the specified pathway names and extracts the corresponding genes based on the selected annotation column.
#'
#' The `annotation` argument allows for flexibility in retrieving genes by different identifiers, such as `gene_symbol` or `gene_id`.
#'
#' @return A vector of genes associated with the specified pathways. The type of genes (e.g., symbols or Ensembl IDs) depends on the `annotation` argument.
#'
#' @importFrom dplyr filter pull
#'
#' @export
#'
get_pathway_genes <- function(annotation_collection,
                              pathway_names,
                              annotation = "gene_symbol") {
  genes <- annotation_collection |>
    dplyr::filter(.data$pathway %in% pathway_names) |>
    dplyr::pull(!! rlang::sym(annotation))

  return(genes)
}
