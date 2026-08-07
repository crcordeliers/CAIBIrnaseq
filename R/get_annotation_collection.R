#' Get Annotation Collection
#'
#' Retrieves gene sets from specified collections for a given species using the MSigDB database.
#'
#' @param collections A character vector specifying the names of the collections to retrieve. Collections can include MSigDB subcategories or "HALLMARK".
#' @param species A character string specifying the species for which the gene sets should be retrieved. Default is `"Homo sapiens"`.
#'
#' @details
#' The function queries the MSigDB database via the `msigdbr` package to collect gene sets for the specified collections and species.
#' If a collection is not found in the MSigDB subcategories, a warning message is displayed, and that collection is skipped.
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
get_annotation_collection <- function(collections, species = "Homo sapiens") {
  # --- INPUT CHECKS ---
  if (missing(collections) || !is.character(collections)) {
    stop("`collections` must be a character vector.")
  }

  if (!is.character(species) || length(species) != 1) {
    stop("`species` must be a single character string.")
  }

  db_species <- case_when(
    species == "Homo sapiens" ~ "HS",
    species == "Mus musculus" ~ "MM",
    TRUE ~ species
  )

  # --- Get available collections from MSigDB ---
  available_collections <- msigdbr::msigdbr_collections(db_species = db_species) |> pull(gs_collection_name, gs_subcollection)
  names(available_collections)[available_collections == "Hallmark"] <- "H"


  collections <- sapply(collections, function(collection) {
    if(collection %in% available_collections) {
      return(collection)
    } else if (collection %in% names(available_collections)) {
      new_collection_name <- available_collections[collection == names(available_collections)]
      return(new_collection_name)
    } else {
      return(NA)
    }
  })

  all_collections <- collections

  results <- lapply(collections, function(collection) {
    if (collection %in% all_collections) {
      message("-- Collecting ", collection, " from MSigDB...")

      msigdb <- msigdbr::msigdbr(species = species, db_species = db_species)

      gene_sets <- msigdb |>
        filter(gs_collection_name %in% collection | gs_collection %in% collection) |>
        mutate(collection = collection) |>
        rename(pathway = gs_name) |>
        select(collection, pathway, gene_id = ensembl_gene, gene_symbol)

      return(gene_sets)
    } else {
      warning("Collection `", collection, "` not found in MSigDB. Skipping.")
      return(NULL)
    }
  })

  # --- Bind and deduplicate ---
  collection_sets <- bind_rows(results)
  if (nrow(collection_sets) == 0) {
    message("No valid gene sets found. Returning NULL.")
    return(NULL)
  }

  return(distinct(collection_sets))
}
