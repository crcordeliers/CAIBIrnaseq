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
#' @param pcutoff adjusted p-value to be used to filter differentially expressed genes.
#' @param direction either "up" or "down". ORA takes only "up" or "down"-regulated
#' genes to do it's overrepresentation analysis.
#'
#' @return A data frame of enriched pathways with columns:
#' \describe{
#'   \item{pathway}{Name of the enriched pathway}
#'   \item{pval}{Raw p-value from Fisher's exact test}
#'   \item{padj}{Adjusted p-value (Benjamini-Hochberg)}
#'   \item{geneRatio}{Proportion of input genes found in the pathway}
#'   \item{bgRatio}{Proportion of background genes found in the pathway}
#'   \item{genes}{Comma-separated list of matched genes}
#' }
#' Column names match the naming used by \code{\link{pathwayFGSEA}} (\code{pathway}, \code{pval}, \code{padj})
#' so that downstream consumers such as \code{\link{plot_pathway_dotplot}} can handle both result types.
#'
#' @importFrom dplyr filter mutate arrange bind_rows desc relocate pull
#' @importFrom tibble rownames_to_column tibble
#' @importFrom stats fisher.test p.adjust
#' @export
#'

pathwayORA <- function(diffexp_result, pathways,
                       direction = "up",
                       id_col = "gene_symbol",
                       pcutoff = 0.05
                       ) {

  message("--- Running Overrepresentation Analysis")

  gene_ids <- rownames(diffexp_result)
  gene_ids <- gene_ids[gene_ids %in% pathways[[id_col]]]

  if (length(gene_ids) == 0) {
    stop("`exp_data` uses unknown gene annotation.
         Try either using ensembl_gene_id or gene_name/gene_symbol.")
  }
  if(direction == "up") {
    genes <- diffexp_result |>
      rownames_to_column("genes") |>
      arrange(padj, desc(log2FoldChange)) |>
      filter(padj < pcutoff, log2FoldChange > 0, genes %in% gene_ids) |>
      pull(genes)
    message("Found ", length(genes), " up-regulated expressed genes...")
  } else if (direction == "down") {
    message("Running ORA on down-regulated genes...")
    genes <- diffexp_result |>
      rownames_to_column("genes") |>
      arrange(padj, log2FoldChange) |>
      filter(padj < pcutoff, log2FoldChange < 0, genes %in% gene_ids) |>
      pull(genes)
    message("Found ", length(genes), " down-regulated expressed genes...")
  }


  # Generate gene set lists
  gene_sets <- split(pathways[[id_col]], pathways$pathway)

  # Universe size
  univ <- length(unique(pathways[[id_col]]))
  gene_ids_filtered <- intersect(gene_ids, unique(pathways[[id_col]]))

  # Perform Over-representation Analysis
  enrich_res <- lapply(names(gene_sets), function(id) {
    path <- gene_sets[[id]]
    genes <- gene_ids_filtered
    ginpath <- length(intersect(genes, path))
    gopath <- length(genes) - ginpath
    opath <- length(path) - ginpath
    rest <- univ - length(path) - gopath


    ctg <- matrix(c(ginpath, opath, gopath, rest), nrow = 2)
    if (any(!is.finite(ctg)) || any(ctg < 0)) {
      stop("The resulting contingency matrix is invalid")
    }

    pfish <- fisher.test(ctg, alternative = "greater")$p.value

    data.frame(
      pathway = id,
      pval = pfish,
      geneRatio = paste0(ginpath, "/", length(genes)),
      bgRatio = paste0(length(path), "/", univ),
      genes = paste(genes[genes %in% path], collapse = ", ")
    )
  }) %>% bind_rows()

  # Adjust p-values and filter
  enrich_res <- enrich_res %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    relocate(padj, .after = pval) %>%
    arrange(padj) %>%
    tibble()

  return(enrich_res)
}
