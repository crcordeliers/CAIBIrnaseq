# Robust Coefficient of Variation

This function computes the robust coefficient of variation for a given
vector. The robust coefficient of variation is calculated as the median
absolute deviation (MAD) divided by the median of the data, providing a
measure of variability that is less sensitive to outliers.

## Usage

``` r
robust_cv(v)
```

## Arguments

- v:

  A numeric vector of values for which the robust coefficient of
  variation is calculated.

## Value

A numeric value representing the robust coefficient of variation for the
input vector.

## Examples

``` r
robust_cv(c(1, 2, 3, 4, 5))
#> [1] 0.4942
robust_cv(c(1, 1, 1, 1, 1))
#> [1] 0
```
