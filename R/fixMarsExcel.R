#' Fix Excel-Corrupted Gene Names
#'
#' This function corrects gene names in a data frame that have been
#' automatically converted by Excel into date formats (e.g., "03/01/01")
#' and replaces them with their original gene names.
#'
#' @param mat A data frame containing a `gene_name` column, possibly with corrupted gene names.
#'
#' @return A data frame with corrected `gene_name` entries.
#' @export
#'
fixMarsExcel <- function(mat) {
  rowMars <- which(mat$gene_name %in% c("03/01/01", "03/01/02"))
  mat$gene_name[rowMars] <- c("MARS1", "MARS2")
  return(mat)
}
