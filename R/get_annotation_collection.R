#' Get Annotation Collection
#'
#' Retrieves gene sets from specified collections for a given species using the MSigDB database.
#'
#' @param collections A character vector specifying the names of the collections to retrieve. Collections can include MSigDB subcategories or "HALLMARKS".
#' @param species A character string specifying the species for which the gene sets should be retrieved. Default is `"Homo sapiens"`.
#'
#' @details
#' The function queries the MSigDB database via the `msigdbr` package to collect gene sets for the specified collections and species.
#' If a collection is not found in the MSigDB subcategories, a warning message is displayed, and that collection is skipped.
#'
#' For the "HALLMARKS" collection, the function maps the "H" category to "HALLMARKS" for consistency.
#'
#' The output contains unique gene sets, including the collection name, pathway name, Ensembl gene IDs, and gene symbols.
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{collection}{The name of the collection.}
#'     \item{pathway}{The pathway name.}
#'     \item{gene_id}{The Ensembl gene ID.}
#'     \item{gene_symbol}{The gene symbol.}
#'   }
#' If no valid collections are found, the function returns `NULL`.
#'
#' @importFrom msigdbr msigdbr msigdbr_collections
#' @importFrom dplyr mutate filter rename select distinct if_else bind_rows
#'
#' @export
#'
get_annotation_collection <- function(collections, species = "Homo sapiens") {
  collection_sets <- lapply(collections, function(collection) {

    if(collection %in% msigdbr::msigdbr_collections()$gs_subcollection | collection == "Hallmark") {
      message("-- Collecting ", collection, " from MSigDB...")

      msigdb <- msigdbr::msigdbr(species = species) |>
        dplyr::mutate(gs_subcollection = dplyr::if_else(gs_collection == "H", "Hallmark", gs_subcollection))

      gene_sets <- msigdb |>
        dplyr::filter(gs_subcollection %in% collection) |>
        dplyr::mutate(collection = collection) |>
        dplyr::rename(pathway = gs_name) |>
        dplyr::select(collection, pathway, gene_id = ensembl_gene, gene_symbol)
      return(gene_sets)
    } else {
      message("Collection `", collection, "` not found.")
      return(NULL)
    }
  }) |> dplyr::bind_rows() |> dplyr::distinct()
  return(collection_sets)
}
