#' Estimate Immune and Stromal Cell Populations using MCPcounter
#'
#' This function estimates the abundance of various immune and stromal cell populations
#' in bulk transcriptomic data using the MCPcounter algorithm (for human) or mMCPcounter (for mouse).
#'
#' @param exp_data A `SummarizedExperiment` object containing normalized gene expression data.
#' @param species Character. Either "Homo sapiens" or "Mus musculus", specifying the organism.
#' @param feature_type Character. Either "gene_name" or "ensembl_gene_id", specifying the gene annotation used in rownames.
#'
#' @return A data frame with estimated cell population scores for each sample.
#' @export
#'
#' @importFrom MCPcounter MCPcounter.estimate
#' @importFrom mMCPcounter mMCPcounter.estimate
#'
mcp_counter <- function(exp_data, species, feature_type = "gene_name") {
  # Extract TPM expression matrix
  gexp <- assays(exp_data)[["tpm"]]

  # Define feature types based on input
  if (feature_type == "gene_name") {
    ft1 <- "HUGO_symbols"
    ft2 <- "Gene.Symbol"
  } else if (feature_type == "ensembl_gene_id") {
    ft1 <- "ENSEMBL_ID"
    ft2 <- "ENSEMBL.ID"
  } else {
    stop("Unsupported feature_type. Choose either 'gene_name' or 'ensembl_gene_id'.")
  }

  # Run MCPcounter for human or mouse
  if (species == "Homo sapiens") {
    mcp_res <- MCPcounter::MCPcounter.estimate(gexp, featuresType = ft1)
  } else if (species == "Mus musculus") {
    mcp_res <- mMCPcounter::mMCPcounter.estimate(gexp, features = ft2)
  } else {
    stop("Unsupported species. Choose either 'Homo sapiens' or 'Mus musculus'.")
  }

  # Return the resulting data frame
  return(mcp_res)
}
