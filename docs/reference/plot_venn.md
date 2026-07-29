# Plot a Venn diagram for two or three vectors

This function plots a Venn diagram using the \`euler\` function from the
eulerr package.

## Usage

``` r
plot_venn(
  v1,
  v2,
  universe_size = NULL,
  v1_name = "V1",
  v2_name = "V2",
  v3 = NULL,
  v3_name = "V3",
  fills = TRUE,
  quantities = TRUE,
  title = NULL,
  ...
)
```

## Arguments

- v1:

  A vector of elements in set 1.

- v2:

  A vector of elements in set 2.

- universe_size:

  Optional total size of the universe for Fisher's test (2-set only).

- v1_name:

  Name for set 1.

- v2_name:

  Name for set 2.

- v3:

  Optional vector of elements in set 3 for a 3-way Venn diagram.

- v3_name:

  Name for set 3.

- fills:

  Colors to fill the diagram.

- quantities:

  Whether to display quantities.

- title:

  Title of the plot.

- ...:

  Additional parameters passed to \`plot()\`.

## Value

A ggplot object representing the Venn diagram.
