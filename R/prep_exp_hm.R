#' Prepare Gene Expression Data for Heatmap Visualization
#'
#' This function extracts and formats gene expression data from a `SummarizedExperiment`
#' object for heatmap plotting. It returns a long-format data frame along with identifiers for samples and features.
#'
#' @param expData A `SummarizedExperiment` object.
#' @param genes Either a character vector of gene names or IDs to extract, or
#' a named list of such vectors (e.g. one per pathway) to produce a
#' row-faceted heatmap, with each list element plotted as its own row group
#' labeled by its name.
#' @param assay Character. The name of the assay to use (default is "norm").
#' @param gene_name Character. The column name in `rowData` containing gene names. Default is "gene_name".
#'
#' @returns A list containing:
#' \describe{
#'   \item{table}{A data frame with expression values and sample annotations.}
#'   \item{colv}{Column variable to use as the x-axis (typically sample IDs).}
#'   \item{rowv}{Row variable(s) corresponding to the genes: a plain vector,
#'   or a named list of vectors when `genes` was a named list.}
#' }
#' @export
#'
#' @importFrom SummarizedExperiment rowData assays colData
#' @importFrom tibble rownames_to_column as_tibble tibble
#' @importFrom dplyr left_join
#'
prep_exp_hm <- function(expData,
                        genes,
                        assay = "norm",
                        gene_name = NULL) {

  genes_list <- if (is.list(genes)) genes else NULL
  genes <- if (is.list(genes)) unique(unlist(genes, use.names = FALSE)) else genes

  gene_annot <- SummarizedExperiment::rowData(expData)
  if(is.null(gene_name)) {
    if (!any(genes %in% rownames(expData))) {
      stop("No `genes` were found in `expData` object. Please check your `genes` input and ensure they match the identifiers in `expData`.")
    }
    missing_genes <- genes[!genes %in% rownames(expData)]
    if (length(missing_genes) > 0) {
      warning(paste0("The following `genes` were not found in `expData` and will be ignored: ", paste(missing_genes, collapse = ", ")))
    genes <- genes[genes %in% rownames(expData)]
    }
  } else {
    if (!gene_name %in% colnames(gene_annot)) {
      stop(paste0("The specified `gene_name` column '", gene_name, "' does not exist in the rowData of `expData`."))
    }
    keep_genes <- gene_annot |>
      as_tibble() |>
      filter(!!sym(gene_name) %in% genes) |>
      pull(!!sym(gene_name))

    missing_genes <- setdiff(genes, keep_genes)
    if (length(missing_genes) > 0) {  
      warning(paste0("The following `genes` were not found in the specified `gene_name` column and will be ignored: ", paste(missing_genes, collapse = ", ")))
    }
  }

  if (!any(genes %in% rownames(expData))) {
    stop("No `genes` were found in `expData` object. Please check your `genes` input and ensure they match the identifiers in `expData`.")
  } else if(!all(genes %in% gene_annot$gene_name)) {
    genes <- genes[genes %in% rownames(expData)]
  } else {
    genes <- genes
  }
  gexp <- assays(expData)[[assay]][genes,]
  gexp_t <- gexp |>
    t() |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    tibble()

  feats <- colnames(gexp_t)[-1]

  samp_annot <- colData(expData) |> as.data.frame()

  gexp_t <- gexp_t |> left_join(samp_annot, by = "sample_id")

  if (is.null(genes_list)) {
    rowv <- feats
  } else {
    # Row-faceted heatmaps require disjoint groups, but gene sets (e.g.
    # pathway leading edges) often overlap; a gene already assigned to an
    # earlier group is dropped from later ones rather than duplicated.
    seen <- character(0)
    rowv <- lapply(genes_list, function(g) {
      g <- setdiff(intersect(g, feats), seen)
      seen <<- c(seen, g)
      g
    })
    dup_dropped <- any(vapply(genes_list, function(g) length(intersect(g, feats)), integer(1)) != lengths(rowv))
    if (dup_dropped) {
      message("Some genes appeared in more than one group; each was kept in the first group it appeared in.")
    }
    empty_groups <- names(rowv)[lengths(rowv) == 0]
    if (length(empty_groups) > 0) {
      warning(paste0("The following gene groups had no matching genes and will be dropped: ", paste(empty_groups, collapse = ", ")))
    }
    rowv <- rowv[lengths(rowv) > 0]
    if (length(rowv) == 0) {
      stop("None of the gene groups in `genes` matched any genes in `expData`.")
    }
  }

  return(list(table = gexp_t, colv = "sample_id", rowv = rowv))
}

