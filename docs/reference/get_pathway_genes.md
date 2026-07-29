# Get Pathway Genes

Extracts genes associated with specific pathways from an annotation
collection.

## Usage

``` r
get_pathway_genes(
  annotation_collection,
  pathway_names,
  annotation = "gene_symbol"
)
```

## Arguments

- annotation_collection:

  A data frame containing pathway annotations. It should have at least
  two columns: \`pathway\` and a gene annotation column (e.g.,
  \`gene_symbol\` or \`gene_id\`).

- pathway_names:

  A character vector specifying the names of pathways for which to
  retrieve associated genes.

- annotation:

  A character string specifying the column in \`annotation_collection\`
  to use for gene annotations. Default is \`"gene_symbol"\`.

## Value

A character vector of genes associated with the specified pathways.

## Details

This function filters the provided annotation collection for the
specified pathway names and extracts the corresponding genes based on
the selected annotation column.
