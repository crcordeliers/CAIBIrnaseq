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
fixMarsExcel <- function(mat) {
  # --- VALIDATION CHECKS ---
  if (!is.data.frame(mat)) {
    stop("Error: `mat` must be a data frame.")
  }

  if (!"gene_name" %in% colnames(mat)) {
    stop("Error: `mat` must contain a `gene_name` column.")
  }

  if (!is.character(mat$gene_name) && !is.factor(mat$gene_name)) {
    stop("Error: `gene_name` column must be of type character or factor.")
  }

  # --- MAIN CORRECTION ---
  corrupted <- c("03/01/01", "03/01/02")
  replacement <- c("MARS1", "MARS2")
  rowMars <- which(mat$gene_name %in% corrupted)

  if (length(rowMars) > 0) {
    mat$gene_name <- as.character(mat$gene_name)  # Ensure it's character
    mat$gene_name[rowMars] <- replacement[match(mat$gene_name[rowMars], corrupted)]
    message("Corrected ", length(rowMars), " Excel-corrupted gene name(s).")
  } else {
    message("No corrupted gene names found.")
  }

  return(mat)
}
