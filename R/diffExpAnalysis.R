#' Differential expression analysis (DESeq2 wrapper)
#'
#' @description
#' Performs a DESeq2-based differential expression workflow on the provided
#' expression data. Automatically constructs the design formula if needed,
#' runs DESeq, extracts results for one or more contrasts, optionally applies
#' LFC shrinkage, and returns sorted tables.
#'
#' @param exp_data A count matrix, SummarizedExperiment containing raw counts
#' and colData for samples.
#' @param design A string or formula (without “~”) specifying the model design
#' (e.g. `"condition"` or `~ condition + batch`). If no leading `~` is found,
#'  it is added automatically.
#' @param lfcShrink Logical; if TRUE (default), performs log₂ fold‑change
#' shrinkage on each result table via DESeq2::lfcShrink.
#' @param contrasts Character vector of named results to extract (as returned
#' by DESeq2::resultsNames). If NULL (default), all names except the intercept
#' are used.
#'
#' @return
#' A named list of data.frames (one per contrast), each sorted by adjusted
#' p‑value then log₂ fold change. If only one contrast is specified, returns
#' a single (unnamed) data.frame.
#'
#' @details
#' - Internally builds or coerces a DESeqDataSet from `exp_data` and `design`.
#' - Runs DESeq2::DESeq on the dataset.
#' - For each contrast, retrieves results with DESeq2::results and optionally
#'   applies shrinkage.
#'
#' @importFrom stringr str_detect
#' @importFrom DESeq2 DESeqDataSet DESeq results lfcShrink resultsNames
#' @export
diffExpAnalysis <- function(exp_data, design, lfcShrink = TRUE, contrasts = NULL) {
  if(!str_detect(design, "~")) {
    design <- paste0("~", design)
  }
  # Apply DESeq model
  ddsSE <- DESeqDataSet(exp_data, design = as.formula(design))
  ddsSE <- DESeq(ddsSE)

  if(is.null(contrasts)) {
    contrasts <- resultsNames(ddsSE)[-1]
  }
  # Get results for all contrasts, applying shrink
  results <- lapply(contrasts, function(contrast) {
    res <- results(ddsSE,
                   name = contrast)
    if(lfcShrink) {
      res <- lfcShrink(ddsSE, coef = contrast, res = res)
    }
    res <- res |>
      as.data.frame() |>
      arrange(padj, desc(log2FoldChange))

    return(res)
  })
  names(results) <- contrasts

  if(length(results) == 1) {
    message("Retuning condition = ", names(results)[1])
    results <- results[[1]]
  }
  return(results)
}

