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
  # Order by decreasing log2FoldChange
  diffexp <- diffexp |>
    dplyr::arrange(desc(log2FoldChange))

  # Named vector of stats
  stat <- as.numeric(diffexp$log2FoldChange)
  names(stat) <- rownames(diffexp)

  # Create list of gene sets from pathway collection
  pathwayList <- split(pathwayCollection$gene_symbol, pathwayCollection$pathway)

  # Run FGSEA
  result <- fgsea::fgseaMultilevel(pathways = pathwayList, stats = stat)

  return(result)
}
