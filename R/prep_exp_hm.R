#' Prepare Gene Expression Data for Heatmap Visualization
#'
#' This function extracts and formats gene expression data from a `SummarizedExperiment`
#' object for heatmap plotting. It returns a long-format data frame along with identifiers for samples and features.
#'
#' @param expData A `SummarizedExperiment` object.
#' @param genes Character vector of gene names or IDs to extract.
#' @param assay Character. The name of the assay to use (default is "norm").
#' @param gene_name Character. The column name in `rowData` containing gene names. Default is "gene_name".
#'
#' @returns A list containing:
#' \describe{
#'   \item{table}{A data frame with expression values and sample annotations.}
#'   \item{colv}{Column variable to use as the x-axis (typically sample IDs).}
#'   \item{rowv}{Row variable(s) corresponding to the genes.}
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

  return(list(table = gexp_t, colv = "sample_id", rowv = feats))
}

