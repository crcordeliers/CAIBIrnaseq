#' @title Pathway Enrichment Using FGSEA
#'
#' @description
#' Performs fast gene set enrichment analysis (FGSEA) using the `fgseaMultilevel` method,
#' based on a ranked list of differentially expressed genes and a collection of pathways.
#'
#' @param diffexp A data frame of differential expression results, with row names as gene IDs
#' and a `log2FoldChange` column used for ranking genes.
#'
#' @param pathwayCollection A data frame with at least two columns: one for pathway names (`pathway`)
#' and one for gene symbols (`gene_symbol`).
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
#' @importFrom dplyr arrange
#' @importFrom fgsea fgseaMultilevel
#' @export
pathwayFGSEA <- function(diffexp, pathwayCollection) {
  # Validation des entrées
  if (!"log2FoldChange" %in% colnames(diffexp)) {
    stop("The 'diffexp' data frame must contain a 'log2FoldChange' column.")
  }

  if (!all(c("pathway", "gene_symbol") %in% colnames(pathwayCollection))) {
    stop("The 'pathwayCollection' data frame must contain 'pathway' and 'gene_symbol' columns.")
  }

  diffexp <- diffexp |>
    dplyr::arrange(desc(log2FoldChange))

  stat <- as.numeric(diffexp$log2FoldChange)
  names(stat) <- rownames(diffexp)

  pathwayList <- split(pathwayCollection$gene_symbol, pathwayCollection$pathway)

  message("Running FGSEA analysis...")

  result <- fgsea::fgseaMultilevel(pathways = pathwayList, stats = stat)

  result <- result |>
    arrange(padj, NES) |>
    as_tibble()

  return(result)
}
