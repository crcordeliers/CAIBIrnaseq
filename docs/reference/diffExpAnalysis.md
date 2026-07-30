# Differential expression analysis (DESeq2 / limma-voom wrapper)

Performs a differential expression workflow on the provided expression
data, using either DESeq2 or limma-voom. Automatically constructs the
design formula if needed, fits the model, extracts results for one or
more contrasts, and returns sorted tables with a consistent set of
columns regardless of \`method\`.

## Usage

``` r
diffExpAnalysis(
  exp_data,
  design,
  method = c("deseq2", "limma-voom")[1],
  lfcShrink = TRUE,
  contrasts = NULL,
  fname = NULL
)
```

## Arguments

- exp_data:

  A count matrix, SummarizedExperiment containing raw counts and colData
  for samples.

- design:

  A string or formula (without “~”) specifying the model design (e.g.
  \`"condition"\` or \`~ condition + batch\`). If no leading \`~\` is
  found, it is added automatically.

- method:

  Character. Either \`"deseq2"\` (default) or \`"limma-voom"\`. Both
  return the same columns (\`log2FoldChange\`, \`pvalue\`, \`padj\`,
  ...) so downstream functions (\`pathwayORA\`, \`pathwayFGSEA\`,
  \`plot_exp_volcano\`, ...) work unmodified regardless of which engine
  was used.

- lfcShrink:

  Logical; if TRUE (default), performs log₂ fold‑change shrinkage on
  each result table via DESeq2::lfcShrink. Only applies when \`method =
  "deseq2"\`; limma-voom already moderates its statistics via
  \`limma::eBayes()\`, so this is ignored (with a message) for \`method
  = "limma-voom"\`.

- contrasts:

  Character vector of named results to extract. Either a single string
  in \`"factor_level_vs_reference"\` form (as returned by
  \`DESeq2::resultsNames()\`), or a 3-element vector \`c(factor, level1,
  level2)\`. Both forms work identically for either \`method\`. If NULL
  (default), all non-intercept contrasts are used.

- fname:

  if provided, will save a \`csv\` file. If multiple contrasts, will use
  the base of fname and append the contrast for the file names.

## Value

A named list of data.frames (one per contrast), each sorted by adjusted
p‑value then log₂ fold change. If only one contrast is specified,
returns a single (unnamed) data.frame.

## Details

\- \`method = "deseq2"\`: builds/coerces a DESeqDataSet from
\`exp_data\` and \`design\`, runs DESeq2::DESeq, and retrieves results
with DESeq2::results (optionally shrunk via DESeq2::lfcShrink). -
\`method = "limma-voom"\`: builds a design matrix from \`exp_data\` and
\`design\`, runs \`limma::voom()\` + \`limma::lmFit()\` on
TMM-normalized counts (via \`edgeR::calcNormFactors()\`), then extracts
each contrast with \`limma::makeContrasts()\` +
\`limma::contrasts.fit()\` + \`limma::eBayes()\` +
\`limma::topTable()\`. Output columns are renamed to match DESeq2's
(\`logFC\` -\> \`log2FoldChange\`, \`P.Value\` -\> \`pvalue\`,
\`adj.P.Val\` -\> \`padj\`).
