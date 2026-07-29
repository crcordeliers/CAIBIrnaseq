# Pre-process a gene expression matrix to have different gene annotation

Pre-process a gene expression matrix to have different gene annotation

## Usage

``` r
gexp_preprocess(
  gexp,
  gene_annotation = NULL,
  keep_annot = "gene_name",
  og_annot = "gene_id",
  keep_stat = sum
)
```

## Arguments

- gexp:

  A gene expression matrix, with samples in columns and genes in rows

- gene_annotation:

  A data.frame of gene annotation, refering to the rows of \`gexp\`. It
  must contain at least two columns: one for the \`keep_annot\` and one
  for the \`og_annot\`

- keep_annot:

  The name of the column in \`gene_annotation\` that contains the gene
  names that you would like to put as rownames of your \`gexp\`, i.e.
  the kept annotation

- og_annot:

  The name of the column in \`gene_annotation\` that contains the gene
  names currently being used by \`gexp\`, i.e. the original annotation.

- keep_stat:

  a function for the statistical measure used for keeping a gene when
  \`keep_annot\` doesn't map uniquely to \`og_annot\`.

## Value

A processed gene expression matrix
