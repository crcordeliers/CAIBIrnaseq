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
#' @importFrom dplyr %>% top_n
#' @importFrom stats sd
#' @importFrom biomaRt useMart getBM
#'
gexp_preprocess <- function(gexp, gene_annotation = NULL , keep_annot = "gene_name",
                            og_annot = "gene_id", keep_stat = sum) {
  if(is.null(gene_annotation)){
    # Connect to the Ensembl database
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") # Replace with "hsapiens_gene_ensembl" for human data

    #-- Retrieve HUGO gene symbols
    gene_symbols <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
                          filters = 'ensembl_gene_id',
                          values = rownames(gene_mat),
                          mart = ensembl)

    #-- Change ensembl rownames with mgi symbols
    rownames(gene_mat) <- gene_symbols$mgi_symbol[match(rownames(gene_mat), gene_symbols$ensembl_gene_id)]
  }
  else{
    #-- Get important annotation
    samples <- colnames(gexp)
    gannot <- dplyr::select(gene_annotation, {{og_annot}}, {{keep_annot}})

    if(as.character(substitute(keep_stat)) %in% c("cv", "sd", "var")) {
      #-- Calculate the keep_stat
      stat <- apply(gexp, 1, keep_stat)
      gexp_int <- cbind(gexp, stat)

      #-- Filter IDs by the new annotation based on max stat per group
      gexp_int <- gexp_int %>%
        as.data.frame() %>%
        rownames_to_column(og_annot) %>%
        left_join(gannot, by = og_annot) %>%
        dplyr::filter(!is.na(!!sym(keep_annot)), !!sym(keep_annot) != "") %>%
        group_by(!!sym(keep_annot)) %>%
        top_n(1, stat) %>%
        dplyr::slice(1) %>%
        ungroup()

      #-- Re-assemble gexp
      gexp_res <- gexp_int %>%
        dplyr::select(-!!og_annot, -stat) %>%
        as.data.frame() %>%
        tibble::column_to_rownames(keep_annot) %>%
        as.matrix()

    } else {
      gexp_res <- gexp %>%
        as.data.frame() %>%
        tibble::rownames_to_column(og_annot) %>%
        left_join(gannot, by = og_annot) %>%
        dplyr::filter(!is.na(!!sym(keep_annot)), !!sym(keep_annot) != "") %>%
        group_by(!!sym(keep_annot)) %>%
        dplyr::select(-!!sym(og_annot)) %>%
        summarise_all(.funs = keep_stat, na.rm = TRUE) %>%
        ungroup() %>%
        as.data.frame() %>%
        column_to_rownames(keep_annot)
    }
  }
  return(gexp_res)
}


# Coefficient of variation function
cv <- function (v) { stats::sd(v, na.rm = TRUE)/mean(v, na.rm = TRUE) }
