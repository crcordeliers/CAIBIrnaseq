#' Extract Genes from Selected Pathways
#'
#' This function extracts gene identifiers from selected pathways within an annotation collection.
#'
#' @param annotation_collection A data frame containing pathway annotations (e.g., from `get_annotation_collection()`).
#' @param pathway_names A character vector of pathway names to extract genes from.
#' @param annotation The column to use for gene identifiers. Default is "gene_symbol".
#'
#' @returns A character vector of gene identifiers present in the specified pathways.
#' @export
#'
#' @importFrom dplyr filter pull
#'
get_pathway_genes <- function(annotation_collection,
                              pathway_names,
                              annotation = "gene_symbol") {
  genes <- annotation_collection |>
    dplyr::filter(.data$pathway %in% pathway_names) |>
    dplyr::pull(!! rlang::sym(annotation))

  return(genes)
}
