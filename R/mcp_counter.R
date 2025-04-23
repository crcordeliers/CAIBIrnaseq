#' MCP-counter wrapper for SummarizedExperiment objects
#'
#' Estimates immune and stromal cell populations using MCPcounter or mMCPcounter.
#'
#' @param exp_data A SummarizedExperiment object containing TPM expression data in the "tpm" assay.
#' @param species Character. Either "Homo sapiens" or "Mus musculus".
#' @param feature_type Character. Either "gene_name" or "ensembl_gene_id".
#'
#' @return A data.frame with MCPcounter scores.
#' @export
#' @importFrom SummarizedExperiment assays
#'
mcp_counter <- function(exp_data, species, feature_type = "gene_name") {
  if (!requireNamespace("MCPcounter", quietly = TRUE)) {
    stop("The 'MCPcounter' package is required but not installed.")
  }

  # Vérification de l'objet SummarizedExperiment
  if (!"SummarizedExperiment" %in% class(exp_data)) {
    stop("exp_data must be a SummarizedExperiment object.")
  }

  gexp <- SummarizedExperiment::assays(exp_data)[["tpm"]]
  if (is.null(gexp)) {
    stop("The 'tpm' assay is missing in the SummarizedExperiment object.")
  }

  if (feature_type == "gene_name") {
    ft1 <- "HUGO_symbols"
    ft2 <- "Gene.Symbol"
  } else if (feature_type == "ensembl_gene_id") {
    ft1 <- "ENSEMBL_ID"
    ft2 <- "ENSEMBL.ID"
  } else {
    stop("Invalid feature_type. Must be either 'gene_name' or 'ensembl_gene_id'.")
  }

  if (species == "Homo sapiens") {
    mcp_res <- MCPcounter::MCPcounter.estimate(gexp, featuresType = ft1)
  } else if (species == "Mus musculus") {
    mcp_res <- mMCPcounter::mMCPcounter.estimate(gexp, features = ft2)
  } else {
    stop("Unsupported species. Must be 'Homo sapiens' or 'Mus musculus'.")
  }

  return(mcp_res)
}
