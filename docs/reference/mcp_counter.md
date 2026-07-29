# MCP-counter wrapper for SummarizedExperiment objects

Estimates immune and stromal cell populations using MCPcounter or
mMCPcounter.

## Usage

``` r
mcp_counter(exp_data, species, feature_type = "gene_name", assay = "tpm")
```

## Arguments

- exp_data:

  A SummarizedExperiment object containing TPM expression data in the
  "tpm" assay.

- species:

  Character. Either "Homo sapiens" or "Mus musculus".

- feature_type:

  Character. Either "gene_name" or "ensembl_gene_id".

## Value

A data.frame with MCPcounter scores.
