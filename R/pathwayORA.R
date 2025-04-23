#' Pathway Over-Representation Analysis (ORA)
#'
#' @description
#' Performs an over-representation analysis (ORA) using Fisher's exact test
#' on a list of differentially expressed genes against a set of predefined pathways.
#'
#' @param diffexp_result A data frame of differential expression results,
#' where row names correspond to gene identifiers.
#' @param pathways A data frame with at least two columns: one indicating
#' pathway names (e.g., 'pathway') and one with gene identifiers.
#' @param id_col Character. The column name in `pathways` that matches gene identifiers in `diffexp_result`.
#' @param pcutoff Numeric. Adjusted p-value threshold to filter significant pathways (default = 0.05).
#'
#' @return A data frame of enriched pathways with columns:
#' \describe{
#'   \item{Pathway}{Name of the enriched pathway}
#'   \item{PValue}{Raw p-value from Fisher's exact test}
#'   \item{PAdj}{Adjusted p-value (Benjamini-Hochberg)}
#'   \item{GeneRatio}{Proportion of input genes found in the pathway}
#'   \item{BgRatio}{Proportion of background genes found in the pathway}
#'   \item{Genes}{Comma-separated list of matched genes}
#' }
#'
#' @importFrom dplyr filter mutate arrange bind_rows
#' @importFrom stats fisher.test p.adjust
#' @export
#'
pathwayORA <- function(diffexp_result, pathways,
                       id_col = "gene_symbol",
                       pcutoff = 0.05) {

  # Extract gene IDs from differential expression result
  gene_ids <- rownames(diffexp_result)

  # Sanity check: ensure overlap between DE genes and pathway genes
  if (!any(gene_ids %in% pathways[[id_col]])) {
    stop("`diffexp_result` uses unknown gene annotation. ",
         "Try either using ensembl_gene_id or gene_name/gene_symbol.")
  }

  # Create list of genes per pathway
  gene_sets <- split(pathways[[id_col]], pathways$pathway)

  # Define the universe: all unique genes in all pathways
  univ_genes <- unique(unlist(gene_sets))
  univ_size <- length(univ_genes)

  # Perform Over-representation Analysis using Fisher's exact test
  enrich_res <- lapply(names(gene_sets), function(id) {
    path_genes <- unique(gene_sets[[id]])
    sig_genes <- gene_ids

    # Overlap between input genes and pathway
    overlap <- sum(sig_genes %in% path_genes)
    if (overlap == 0) {
      return(NULL)  # Skip pathways with no overlap
    }

    # Construct contingency table
    a <- overlap  # in path and in sig
    b <- length(path_genes) - a  # in path and not in sig
    c <- length(sig_genes) - a  # not in path and in sig
    d <- univ_size - a - b - c  # not in path and not in sig

    # Avoid errors in Fisher test due to 0-counts or invalid counts
    if (any(c(a, b, c, d) < 0) || any(is.infinite(c(a, b, c, d)))) {
      return(NULL)  # Skip this pathway if the contingency table is invalid
    }

    # Construct contingency table
    ctg <- matrix(c(a, b, c, d), nrow = 2)
    pfish <- tryCatch({
      fisher.test(ctg, alternative = "greater")$p.value
    }, error = function(e) {
      return(NA)  # Return NA if fisher.test fails
    })

    data.frame(
      Pathway = id,
      PValue = pfish,
      GeneRatio = paste0(a, "/", length(sig_genes)),
      BgRatio = paste0(b, "/", univ_size - length(path_genes)),
      Genes = paste(sig_genes[sig_genes %in% path_genes], collapse = ", ")
    )
  }) %>% dplyr::bind_rows()

  # Adjust p-values, sort, and filter by cutoff
  if (nrow(enrich_res) > 0) {
    enrich_res <- enrich_res %>%
      dplyr::mutate(PAdj = p.adjust(PValue, method = "BH")) %>%
      dplyr::arrange(PAdj) %>%
      dplyr::filter(PAdj < pcutoff)
  } else {
    enrich_res <- data.frame()  # return empty data frame if no results
  }

  return(enrich_res)
}
