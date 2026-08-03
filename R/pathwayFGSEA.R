#' @title Pathway Enrichment Using FGSEA
#'
#' @description
#' Performs fast gene set enrichment analysis (FGSEA) using the `fgseaMultilevel` method,
#' based on a ranked list of differentially expressed genes and a collection of pathways.
#'
#' @param diffexp A data frame of differential expression results, with row names as gene IDs
#' and a `log2FoldChange` column used for ranking genes.
#'
#' @param pathwayCollection Either a data frame with at least two columns, one for pathway names
#' (`pathway`) and one for gene symbols (`gene_symbol`), or a named list of character vectors
#' (one gene-symbol vector per pathway/signature), e.g. `list(MySignature = c("GENE1", "GENE2"))`.
#' The list form lets custom gene signatures be tested the same way as an MSigDB collection.
#'
#' @return A data frame returned by `fgseaMultilevel()`, including columns such as:
#' \describe{
#'   \item{pathway}{Name of the pathway}
#'   \item{pval}{P-value of enrichment}
#'   \item{padj}{Adjusted p-value (FDR)}
#'   \item{ES}{Enrichment Score}
#'   \item{NES}{Normalized Enrichment Score}
#'   \item{leadingEdge}{Vector of leading-edge genes}
#' }
#'
#' @importFrom dplyr arrange desc mutate pull filter
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom fgsea fgseaMultilevel
#' @export
pathwayFGSEA <- function(diffexp, pathwayCollection, seed = 0) {
  # Validate inputs
  if (!"log2FoldChange" %in% colnames(diffexp)) {
    stop("The 'diffexp' data frame must contain a 'log2FoldChange' column.")
  }

  pathwayList <- .as_gene_sets(pathwayCollection, id_col = "gene_symbol")

  if(!is.null(seed)) {
    set.seed(seed)
  }
  # Genes with NA pvalue/log2FoldChange (e.g. filtered out internally by DESeq2's
  # independent filtering or Cook's-distance outlier detection) can't be ranked
  # and are dropped rather than passed through as non-finite values.
  pheno <- diffexp |>
    rownames_to_column("symbol") |>
    filter(!is.na(pvalue), !is.na(log2FoldChange)) |>
    mutate(stat = -log10(pvalue)*sign(log2FoldChange)) |>
    filter(is.finite(stat)) |>
    arrange(stat, desc(stat)) |>
    pull(stat, symbol)

  message("Running FGSEA analysis...")

  result <- fgseaMultilevel(pathways = pathwayList, pheno)

  result <- result |>
    arrange(padj, NES) |>
    as_tibble()

  return(result)
}
