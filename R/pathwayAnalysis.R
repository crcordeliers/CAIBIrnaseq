#' Pathway Analysis
#'
#' Performs pathway enrichment analysis using either Over-representation Analysis (ORA)
#' or Gene Set Enrichment Analysis (FGSEA) based on the provided gene set and pathways.
#'
#' @param gene_set A data frame containing the differential expression results. Rows represent genes.
#' @param pathways A data frame containing pathway information, with columns for pathways and associated genes.
#' @param method A string specifying the pathway analysis method. Options are `"ORA"` for Over-representation Analysis or `"FGSEA"` for Gene Set Enrichment Analysis. Default is `"ORA"`.
#' @param species A string specifying the species for which the pathways are relevant. Default is `"Homo sapiens"`.
#' @param id_col A string specifying the column in the `pathways` data frame that contains gene identifiers. Default is `"gene_symbol"`.
#' @param pcutoff A numeric value specifying the adjusted p-value cutoff for significant pathways. Default is `0.05`.
#' @param verbose A logical value indicating whether to display progress messages. Default is `TRUE`.
#'
#' @details
#' This function allows users to perform pathway enrichment analysis using either ORA or FGSEA.
#' The `pathways` input should include a column specifying gene identifiers and a column specifying pathway names.
#'
#' - **ORA**: Identifies pathways over-represented in the provided gene set.
#' - **FGSEA**: Identifies pathways using a ranked list of genes based on differential expression statistics.
#'
#' The function dynamically calls internal implementations of ORA or FGSEA based on the `method` argument.
#'
#' @return A data frame containing pathway enrichment results, including pathway names, p-values, adjusted p-values, and other relevant statistics.
#'
#' @export
pathwayAnalysis <- function(gene_set, pathways,
                            method = "ORA",
                            species = "Homo sapiens",
                            id_col = "gene_symbol",
                            pcutoff = 0.05,
                            verbose = TRUE) {

  if (nrow(pathways) == 0) {
    stop("No pathways retrieved from the provided collections.")
  }

  # Select analysis method
  if (tolower(method) == "ora") {
    result <- pathwayORA(gene_set, pathways, id_col = id_col, pcutoff = pcutoff)
  } else if (tolower(method) == "fgsea") {
    result <- pathwayFGSEA(gene_set, pathways)
  } else {
    stop("Invalid method. Use a method in 'ORA' or 'FGSEA'.")
  }
  return(result)
}
