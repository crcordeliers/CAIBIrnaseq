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

  gene_ids <- rownames(diffexp_result)

  if (!any(gene_ids %in% pathways[[id_col]])) {
    stop("`exp_data` uses unknown gene annotation.
         Try either using ensembl_gene_id or gene_name/gene_symbol.")
  }

  # Generate gene set lists
  gene_sets <- split(pathways[[id_col]], pathways$pathway)

  # Universe size
  univ <- length(unique(unlist(gene_sets)))

  # Perform Over-representation Analysis
  enrich_res <- lapply(names(gene_sets), function(id) {
    path <- gene_sets[[id]]
    genes <- rownames(diffexp_result)
    ginpath <- sum(genes %in% path)
    gopath <- length(genes) - ginpath
    opath <- length(path) - ginpath
    rest <- univ - ginpath - gopath - opath
    ctg <- matrix(c(ginpath, opath, gopath, rest), nrow = 2)
    pfish <- fisher.test(ctg, alternative = "greater")$p.value

    data.frame(
      Pathway = id,
      PValue = pfish,
      GeneRatio = paste0(ginpath, '/', length(genes)),
      BgRatio = paste0(opath, '/', univ - length(path)),
      Genes = paste(genes[genes %in% path], collapse = ", ")
    )
  }) %>% bind_rows()

  # Adjust p-values and filter
  enrich_res <- enrich_res %>%
    mutate(PAdj = p.adjust(PValue, method = "BH")) %>%
    arrange(PAdj) %>%
    filter(PAdj < pcutoff)

  return(enrich_res)
}
