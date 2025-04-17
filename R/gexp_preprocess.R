require(tidyverse)

#' Pre-process a gene expression matrix to have different gene annotation
#'
#' @param gexp A gene expression matrix, with samples in columns and genes in rows
#' @param gene_annotation A data.frame of gene annotation, refering to the rows of
#' `gexp`. It must contain at least two columns: one for the `keep_annot` and one
#' for the `og_annot`
#' @param keep_annot The name of the column in `gene_annotation` that contains
#' the gene names that you would like to put as rownames of your `gexp`, i.e. the
#' kept annotation
#' @param og_annot The name of the column in `gene_annotation` that contains the
#' gene names currently being used by `gexp`, i.e. the original annotation.
#' @param keep_stat a function for the statistical measure
#' used for keeping a gene when `keep_annot` doesn't map uniquely to `og_annot`.
#' @returns A processed gene expression matrix
#' @export
#'
#' @importFrom dplyr select filter group_by slice_max ungroup summarise_all across left_join slice everything
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang sym
#' @importFrom dplyr %>%
#' @importFrom stats sd
#'
gexp_preprocess <- function(gexp, gene_annotation, keep_annot = "gene_name",
                            og_annot = "gene_id", keep_stat = sum) {

  samples <- colnames(gexp)
  gannot <- dplyr::select(gene_annotation, {{og_annot}}, {{keep_annot}})

  if (as.character(substitute(keep_stat)) %in% c("cv", "sd", "var")) {
    stat <- apply(gexp, 1, keep_stat)
    gexp_int <- cbind(gexp, stat) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(og_annot) %>%
      dplyr::left_join(gannot, by = og_annot) %>%
      dplyr::filter(!is.na(!!sym(keep_annot)), !!sym(keep_annot) != "") %>%
      dplyr::group_by(!!sym(keep_annot)) %>%
      dplyr::slice_max(order_by = stat, n = 1) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    gexp_res <- gexp_int %>%
      dplyr::select(-!!sym(og_annot), -stat) %>%
      tibble::column_to_rownames(keep_annot) %>%
      as.matrix()

  } else {
    gexp_res <- gexp %>%
      as.data.frame() %>%
      tibble::rownames_to_column(og_annot) %>%
      dplyr::left_join(gannot, by = og_annot) %>%
      dplyr::filter(!is.na(!!sym(keep_annot)), !!sym(keep_annot) != "") %>%
      group_by(!!sym(keep_annot)) %>%
      dplyr::select(-!!sym(og_annot)) %>%
      dplyr::summarise_all(.funs = keep_stat, na.rm = TRUE) %>%
      ungroup() %>%
      as.data.frame() %>%
      tibble::column_to_rownames(keep_annot)
  }

  return(gexp_res)
}


# Coefficient of variation function
cv <- function (v) { stats::sd(v, na.rm = TRUE)/mean(v, na.rm = TRUE) }
