# Get Annotation Collection

Retrieves gene sets from specified collections for a given species using
the MSigDB database.

## Usage

``` r
get_annotation_collection(collections, species = "Homo sapiens")
```

## Arguments

- collections:

  A character vector specifying the names of the collections to
  retrieve. Collections can include MSigDB subcategories or "HALLMARK".

- species:

  A character string specifying the species for which the gene sets
  should be retrieved. Default is \`"Homo sapiens"\`.

## Value

A data frame with the following columns:

- collection:

  The name of the collection.

- pathway:

  The pathway name.

- gene_id:

  The Ensembl gene ID.

- gene_symbol:

  The gene symbol.

If no valid collections are found, the function returns \`NULL\`.

## Details

The function queries the MSigDB database via the \`msigdbr\` package to
collect gene sets for the specified collections and species. If a
collection is not found in the MSigDB subcategories, a warning message
is displayed, and that collection is skipped.
