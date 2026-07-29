# Differential expression analysis (DESeq2 wrapper)

Performs a DESeq2-based differential expression workflow on the provided
expression data. Automatically constructs the design formula if needed,
runs DESeq, extracts results for one or more contrasts, optionally
applies LFC shrinkage, and returns sorted tables.

## Usage

``` r
diffExpAnalysis(
  exp_data,
  design,
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

- lfcShrink:

  Logical; if TRUE (default), performs log₂ fold‑change shrinkage on
  each result table via DESeq2::lfcShrink.

- contrasts:

  Character vector of named results to extract (as returned by
  DESeq2::resultsNames). If NULL (default), all names except the
  intercept are used.

- fname:

  if provided, will save a \`csv\` file. If multiple contrasts, will use
  the base of fname and append the contrast for the file names.

## Value

A named list of data.frames (one per contrast), each sorted by adjusted
p‑value then log₂ fold change. If only one contrast is specified,
returns a single (unnamed) data.frame.

## Details

\- Internally builds or coerces a DESeqDataSet from \`exp_data\` and
\`design\`. - Runs DESeq2::DESeq on the dataset. - For each contrast,
retrieves results with DESeq2::results and optionally applies shrinkage.
