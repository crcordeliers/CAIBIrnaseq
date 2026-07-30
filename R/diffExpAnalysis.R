#' Differential expression analysis (DESeq2 / limma-voom wrapper)
#'
#' @description
#' Performs a differential expression workflow on the provided expression data,
#' using either DESeq2 or limma-voom. Automatically constructs the design
#' formula if needed, fits the model, extracts results for one or more
#' contrasts, and returns sorted tables with a consistent set of columns
#' regardless of `method`.
#'
#' @param exp_data A count matrix, SummarizedExperiment containing raw counts
#' and colData for samples.
#' @param design A string or formula (without “~”) specifying the model design
#' (e.g. `"condition"` or `~ condition + batch`). If no leading `~` is found,
#'  it is added automatically.
#' @param method Character. Either `"deseq2"` (default) or `"limma-voom"`.
#' Both return the same columns (`log2FoldChange`, `pvalue`, `padj`, ...) so
#' downstream functions (`pathwayORA`, `pathwayFGSEA`, `plot_exp_volcano`, ...)
#' work unmodified regardless of which engine was used.
#' @param lfcShrink Logical; if TRUE (default), performs log₂ fold‑change
#' shrinkage on each result table via DESeq2::lfcShrink. Only applies when
#' `method = "deseq2"`; limma-voom already moderates its statistics via
#' `limma::eBayes()`, so this is ignored (with a message) for `method = "limma-voom"`.
#' @param contrasts Character vector of named results to extract. Either a
#' single string in `"factor_level_vs_reference"` form (as returned by
#' `DESeq2::resultsNames()`), or a 3-element vector `c(factor, level1, level2)`.
#' Both forms work identically for either `method`. If NULL (default), all
#' non-intercept contrasts are used.
#' @param fname if provided, will save a `csv` file. If multiple contrasts,
#' will use the base of fname and append the contrast for the file names.
#'
#' @return
#' A named list of data.frames (one per contrast), each sorted by adjusted
#' p‑value then log₂ fold change. If only one contrast is specified, returns
#' a single (unnamed) data.frame.
#'
#' @details
#' - `method = "deseq2"`: builds/coerces a DESeqDataSet from `exp_data` and
#'   `design`, runs DESeq2::DESeq, and retrieves results with DESeq2::results
#'   (optionally shrunk via DESeq2::lfcShrink).
#' - `method = "limma-voom"`: builds a design matrix from `exp_data` and
#'   `design`, runs `limma::voom()` + `limma::lmFit()` on TMM-normalized counts
#'   (via `edgeR::calcNormFactors()`), then extracts each contrast with
#'   `limma::makeContrasts()` + `limma::contrasts.fit()` + `limma::eBayes()` +
#'   `limma::topTable()`. Output columns are renamed to match DESeq2's
#'   (`logFC` -> `log2FoldChange`, `P.Value` -> `pvalue`, `adj.P.Val` -> `padj`).
#'
#' @importFrom stringr str_detect
#' @importFrom DESeq2 DESeqDataSet DESeq results lfcShrink resultsNames
#' @importFrom readr write_csv
#' @importFrom fs path_ext_remove
#' @importFrom dplyr arrange desc rename
#' @importFrom stats model.matrix as.formula
#' @export
diffExpAnalysis <- function(exp_data, design, method = c("deseq2", "limma-voom")[1],
                             lfcShrink = TRUE, contrasts = NULL,
                             fname = NULL) {
  method <- tolower(method)
  if (!method %in% c("deseq2", "limma-voom")) {
    stop("`method` must be either \"deseq2\" or \"limma-voom\".")
  }

  if(!str_detect(design, "~")) {
    design <- paste0("~", design)
  }
  # Coerce counts to integer (required by DESeq2; harmless for limma-voom)
  counts <- SummarizedExperiment::assay(exp_data, "counts")
  if (!is.integer(counts)) {
    counts <- matrix(
      as.integer(round(counts)), nrow = nrow(counts), dimnames = dimnames(counts)
    )
    SummarizedExperiment::assay(exp_data, "counts") <- counts
  }

  if (method == "deseq2") {
    results_list <- .diffExp_deseq2(exp_data, design, lfcShrink, contrasts)
  } else {
    if (isTRUE(lfcShrink)) {
      message("`lfcShrink` is ignored for method = \"limma-voom\" (limma::eBayes() already moderates statistics).")
    }
    results_list <- .diffExp_limma_voom(exp_data, design, contrasts)
  }

  results <- results_list$results
  contrast_names <- results_list$contrast_names
  names(results) <- contrast_names

  if(length(results) == 1) {
    results <- results[[1]]
    if(!is.null(fname)) {
      res <- results |> rownames_to_column("gene")
      write_csv(res, fname)
    }
  } else {
    if(!is.null(fname)) {
      fnames <- sapply(contrast_names, function(ctr) {
        paste0(path_ext_remove(fname), "_", ctr, ".csv") })
      lapply(1:length(results), function(i) {
        res <- results[[i]] |> rownames_to_column("gene")
        write_csv(res, file =  fnames[i])
      })
    }
  }
  return(results)
}

# --- DESeq2 engine ----------------------------------------------------------

.diffExp_deseq2 <- function(exp_data, design, lfcShrink, contrasts) {
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
  contrast_names <- lapply(contrasts, paste, collapse = "_") |> unlist()

  list(results = results, contrast_names = contrast_names)
}

# --- limma-voom engine -------------------------------------------------------

.diffExp_limma_voom <- function(exp_data, design, contrasts) {
  colData_df <- as.data.frame(SummarizedExperiment::colData(exp_data))
  factor_names <- intersect(all.vars(as.formula(design)), colnames(colData_df))

  counts <- SummarizedExperiment::assay(exp_data, "counts")
  design_matrix <- model.matrix(as.formula(design), data = colData_df)

  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge, design_matrix)
  fit <- limma::lmFit(v, design_matrix)

  if(is.null(contrasts)) {
    contrast_triples <- lapply(colnames(design_matrix)[-1], function(coef) {
      .coef_to_triple(coef, factor_names, colData_df)
    })
  } else {
    # Mirrors the DESeq2 branch's semantics: lapply() over `contrasts` treats
    # each *element* as one contrast spec. A flat length-3 character vector is
    # therefore 3 separate single-name contrasts; to request one explicit
    # c(factor, level1, level2) contrast, wrap it in a list.
    contrast_triples <- lapply(contrasts, .as_contrast_triple, factor_names = factor_names, colData_df = colData_df)
  }

  results <- lapply(contrast_triples, function(triple) {
    factor_name <- triple[1]; level1 <- triple[2]; level2 <- triple[3]

    if(identical(level1, "") && identical(level2, "")) {
      # Non-categorical (e.g. numeric) covariate: extract its single coefficient directly.
      fit2 <- limma::eBayes(fit)
      res <- limma::topTable(fit2, coef = factor_name, number = Inf, sort.by = "none")
    } else {
      ref_level <- levels(factor(colData_df[[factor_name]]))[1]
      coef1 <- if(level1 == ref_level) "0" else paste0(factor_name, level1)
      coef2 <- if(level2 == ref_level) "0" else paste0(factor_name, level2)
      missing_coefs <- setdiff(c(coef1, coef2), c("0", colnames(design_matrix)))
      if(length(missing_coefs) > 0) {
        stop("Contrast level(s) not found in the design matrix: ", paste(missing_coefs, collapse = ", "),
             ". Available coefficients: ", paste(colnames(design_matrix), collapse = ", "))
      }

      # `makeContrasts()` warns about renaming the "(Intercept)" design-matrix column to a
      # syntactically valid name; that's an internal housekeeping detail, not user-actionable.
      cm <- suppressWarnings(
        limma::makeContrasts(contrasts = paste(coef1, "-", coef2), levels = design_matrix)
      )
      fit2 <- limma::contrasts.fit(fit, cm)
      fit2 <- limma::eBayes(fit2)
      res <- limma::topTable(fit2, number = Inf, sort.by = "none")
    }

    res <- res |>
      rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val) |>
      arrange(padj, desc(log2FoldChange))

    return(res)
  })

  contrast_names <- vapply(contrast_triples, function(triple) {
    if(identical(triple[2], "") && identical(triple[3], "")) {
      triple[1]
    } else {
      paste(triple[1], triple[2], "vs", triple[3], sep = "_")
    }
  }, character(1))

  list(results = results, contrast_names = contrast_names)
}

# Parse a DESeq2-style "factor_level1_vs_level2" string (or pass through a
# length-3 c(factor, level1, level2) vector) into a c(factor, level1, level2)
# triple. Uses the known design factor names to resolve factor names that
# themselves contain underscores (e.g. "condition_name").
.as_contrast_triple <- function(contrast, factor_names, colData_df) {
  if(length(contrast) == 3) {
    return(as.character(contrast))
  }
  if(length(contrast) != 1) {
    stop("Can't work with the `contrasts` provided.")
  }
  candidates <- factor_names[startsWith(contrast, paste0(factor_names, "_"))]
  if(length(candidates) == 0) {
    stop("Could not match contrast \"", contrast, "\" to a design factor. ",
         "Use a length-3 vector c(factor, level1, level2) instead.")
  }
  factor_name <- candidates[which.max(nchar(candidates))]
  rest <- sub(paste0("^", factor_name, "_"), "", contrast)
  parts <- strsplit(rest, "_vs_")[[1]]
  if(length(parts) != 2) {
    stop("Could not parse levels from contrast \"", contrast, "\". ",
         "Expected \"", factor_name, "_<level1>_vs_<level2>\".")
  }
  c(factor_name, parts[1], parts[2])
}

# Reverse-map a limma design-matrix coefficient name (e.g. "dextrt") back to
# a c(factor, level, reference) triple, for auto-derived (contrasts = NULL) contrasts.
.coef_to_triple <- function(coef, factor_names, colData_df) {
  candidates <- factor_names[startsWith(coef, factor_names)]
  if(length(candidates) == 0) {
    return(c(coef, "", ""))
  }
  factor_name <- candidates[which.max(nchar(candidates))]
  level <- sub(paste0("^", factor_name), "", coef)
  ref_level <- levels(factor(colData_df[[factor_name]]))[1]
  c(factor_name, level, ref_level)
}
