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
#' @param fname if provided, will save a `csv` file. If multiple contrasts,
#' will use the base of fname and append the contrast for the file names.
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
#' @importFrom readr write_csv
#' @importFrom fs path_ext_remove
#' @export
diffExpAnalysis <- function(exp_data, design, lfcShrink = TRUE, contrasts = NULL,
                            fname = NULL) {
  if(!str_detect(design, "~")) {
    design <- paste0("~", design)
  }
  # Apply DESeq model
  counts <- SummarizedExperiment::assay(exp_data, "counts")
  if (!is.integer(counts)) {
    SummarizedExperiment::assay(exp_data, "counts") <- matrix(
      as.integer(round(counts)), nrow = nrow(counts), dimnames = dimnames(counts)
    )
  }
  ddsSE <- DESeqDataSet(exp_data, design = as.formula(design))
  ddsSE <- DESeq(ddsSE)

  if(is.null(contrasts)) {
    contrasts <- resultsNames(ddsSE)[-1]
  }

  results <- lapply(contrasts, function(contrast) {
    if(length(contrast) == 3) {
      res <- results(ddsSE, contrast = contrast)
      if(lfcShrink) {
        res <- lfcShrink(ddsSE, contrast = contrast, type = "normal", res = res)
      }
    } else if(length(contrast) == 1) {
      res <- results(ddsSE, name = contrast)
      if(lfcShrink) {
        res <- lfcShrink(ddsSE, coef = contrast, res = res)
      }
    } else {
      stop("Can't work with the `contrasts` provided.")
    }
    res <- res |>
      as.data.frame() |>
      arrange(padj, desc(log2FoldChange))

    return(res)
  })
  contrasts <- lapply(contrasts, paste, collapse = "_") |> unlist()

  names(results) <- contrasts

  if(length(results) == 1) {
    results <- results[[1]]
    if(!is.null(fname)) {
      res <- results |> rownames_to_column("gene")
      write_csv(res, fname)
    }
  } else {
    if(!is.null(fname)) {
      fnames <- sapply(contrasts, function(ctr) {
        paste0(path_ext_remove(fname), "_", ctr, ".csv") })
      lapply(1:length(results), function(i) {
        res <- results[[i]] |> rownames_to_column("gene")
        write_csv(res, file =  fnames[i])
      })
    }
  }
  return(results)
}

