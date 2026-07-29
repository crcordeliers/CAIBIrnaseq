#' Plot a Venn diagram for two or three vectors
#'
#' This function plots a Venn diagram using the `euler` function from the eulerr package.
#'
#' @param v1 A vector of elements in set 1.
#' @param v2 A vector of elements in set 2.
#' @param universe_size Optional total size of the universe for Fisher's test (2-set only).
#' @param v1_name Name for set 1.
#' @param v2_name Name for set 2.
#' @param v3 Optional vector of elements in set 3 for a 3-way Venn diagram.
#' @param v3_name Name for set 3.
#' @param fills Colors to fill the diagram.
#' @param quantities Whether to display quantities.
#' @param title Title of the plot.
#' @param ... Additional parameters passed to `plot()`.
#'
#' @return A ggplot object representing the Venn diagram.
#' @importFrom eulerr euler
#' @importFrom ggplotify as.ggplot
#' @importFrom ggplot2 labs
#' @importFrom stats fisher.test
#' @export
plot_venn <- function(v1, v2, universe_size = NULL,
                      v1_name = 'V1', v2_name = 'V2',
                      v3 = NULL, v3_name = 'V3',
                      fills = TRUE, quantities = TRUE,
                      title = NULL, ...) {

  if (!is.null(v3)) {
    # 3-way Venn
    only1   <- length(setdiff(v1, union(v2, v3)))
    only2   <- length(setdiff(v2, union(v1, v3)))
    only3   <- length(setdiff(v3, union(v1, v2)))
    i12     <- length(setdiff(intersect(v1, v2), v3))
    i13     <- length(setdiff(intersect(v1, v3), v2))
    i23     <- length(setdiff(intersect(v2, v3), v1))
    i123    <- length(intersect(intersect(v1, v2), v3))

    n12 <- paste0(v1_name, "&", v2_name)
    n13 <- paste0(v1_name, "&", v3_name)
    n23 <- paste0(v2_name, "&", v3_name)
    n123 <- paste0(v1_name, "&", v2_name, "&", v3_name)

    info4euler <- structure(
      c(only1, only2, only3, i12, i13, i23, i123),
      names = c(v1_name, v2_name, v3_name, n12, n13, n23, n123)
    )

    subtitle <- NULL
    caption  <- paste0("Common to all three: ", i123, " genes")

  } else {
    # 2-way Venn
    v1_only <- length(setdiff(v1, v2))
    v2_only <- length(setdiff(v2, v1))
    v1v2    <- length(intersect(v1, v2))
    v1v2_name <- paste0(v1_name, "&", v2_name)

    info4euler <- structure(c(v1_only, v2_only, v1v2),
                            names = c(v1_name, v2_name, v1v2_name))

    subtitle <- NULL
    if (!is.null(universe_size)) {
      ctg_matrix <- matrix(c(v1v2, v1_only, v2_only, universe_size - sum(v1_only, v2_only, v1v2)),
                           ncol = 2, byrow = TRUE)
      pval <- stats::fisher.test(ctg_matrix)$p.value
      subtitle <- paste0("Fisher test p-value = ", signif(pval, 3))
    }
    p1   <- signif(v1_only / universe_size, 3)
    p1_2 <- signif(v1v2 / universe_size, 3)
    p2   <- signif(v2_only / universe_size, 3)

    caption <- paste0("Proportion of DE significant ", v1_name, " genes: ", p1,
                      "\nProportion of DE significant ", v1v2_name, " genes: ", p1_2,
                      "\nProportion of DE significant ", v2_name, " genes: ", p2)
  }

  plt <- invisible(plot(eulerr::euler(info4euler), fills = fills, quantities = quantities, ...))
  plt <- ggplotify::as.ggplot(plt)
  plt + ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
}

