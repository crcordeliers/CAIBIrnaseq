# Normalize Counts to Reads Per Million (RPM)

This function normalizes raw read counts to Reads Per Million (RPM)
across each column (typically samples). It is commonly used for
comparing sequencing depth-normalized gene expression across samples.

## Usage

``` r
readsPerMillion(data, factor = 10^6)
```

## Arguments

- data:

  A numeric matrix or data frame of raw counts. Rows are typically
  genes, and columns are samples.

- factor:

  A numeric value for scaling. Default is `10^6` to compute reads per
  million.

## Value

A matrix of the same dimensions as `data`, with RPM-normalized values.
