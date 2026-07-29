# Fix Excel-Corrupted Gene Names

This function corrects gene names in a data frame that have been
automatically converted by Excel into date formats (e.g., "03/01/01")
and replaces them with their original gene names.

## Usage

``` r
fixMarsExcel(mat)
```

## Arguments

- mat:

  A data frame containing a \`gene_name\` column, possibly with
  corrupted gene names.

## Value

A data frame with corrected \`gene_name\` entries.
