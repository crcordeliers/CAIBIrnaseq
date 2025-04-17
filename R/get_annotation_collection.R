#' Get Annotation Collections from MSigDB
#'
#' This function retrieves gene sets from specified MSigDB collections.
#'
#' @param collections A character vector of MSigDB sub-collections (e.g., "Hallmark", "CP:KEGG", etc.).
#' @param species The species for which to retrieve gene sets (default: "Homo sapiens").
#'
#' @return A data frame with gene set annotations, including collection name, pathway, Ensembl gene ID, and gene symbol.
#' @export
#'
#' @importFrom msigdbr msigdbr msigdbr_collections
#' @importFrom dplyr mutate if_else filter rename select bind_rows distinct
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
