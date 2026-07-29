# Rebase Gene Expression Matrix

This function rebases the gene expression matrix in a
\`SummarizedExperiment\` object using a specified annotation as the
primary identifier.

## Usage

``` r
rebase_gexp(
  exp_data,
  gene_id = "gene_id",
  new_gene_id = "gene_name",
  keep_cols = NULL
)
```

## Arguments

- exp_data:

  A \`SummarizedExperiment\` object containing the gene expression data.

- gene_id:

  A character string specifying the column in the gene annotation with
  the old main annotation. Default is \`"gene_id"\`.

- new_gene_id:

  A character string specifying the column in the gene annotation to use
  as the main identifier. Default is \`"gene_name"\`.

- keep_cols:

  A character vector specifying columns to keep in the new annotation.
  If \`NULL\`, all original columns are kept and aggregated as
  follows: - For numeric variables, the mean is kept for multimapping
  resolution - For character variables, all values are kept in a long
  string

## Value

A \`SummarizedExperiment\` object with rebased gene expression data. The
output includes updated \`rowData\` with summarized gene-level
information, and new \`assays\` containing rebased \`counts\` and
\`tpm\` matrices.

## Details

This function rebases the gene expression data by aggregating counts
based on the specified annotation. It also calculates transcripts per
million (TPM) using the aggregated data and associated gene lengths.

The function performs the following steps: - Aggregates gene counts and
metadata based on the specified annotation. - If \`gene_lengths_kb\` is
available, calculates TPM values using the rebased gene counts and
average gene lengths. - Constructs a new \`SummarizedExperiment\` object
with the rebased data.
