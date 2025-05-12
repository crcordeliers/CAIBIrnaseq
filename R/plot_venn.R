#' Plot a Venn diagram for two vectors
#'
#' This function plots a Venn diagram using the `euler` function from the eulerr package.
#'
#' @param v1 A vector of elements in set 1.
#' @param v2 A vector of elements in set 2.
#' @param universe_size Optional total size of the universe for Fisher's test.
#' @param v1_name Name for set 1.
#' @param v2_name Name for set 2.
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
                      fills = TRUE, quantities = TRUE,
                      title = NULL, ...) {

  v1_only <- length(setdiff(v1, v2))
  v2_only <- length(setdiff(v2, v1))
  v1v2 <- length(intersect(v1, v2))

  v1v2_name <- paste0(v1_name, "&", v2_name)

  info4euler <- structure(c(v1_only, v2_only, v1v2),
                          names = c(v1_name, v2_name, v1v2_name))

  subtitle <- NULL
  if(!is.null(universe_size)) {
    ctg_matrix <- matrix(c(v1v2, v1_only, v2_only, universe_size - sum(v1_only, v2_only, v1v2)),
                         ncol = 2, byrow = TRUE)
    pval <- stats::fisher.test(ctg_matrix)$p.value
    subtitle <- paste0("Fisher test p-value = ", signif(pval, 3))
  }
  p1 <- signif(v1_only/universe_size, 3)
  p1_2 <- signif(v1v2/universe_size, 3)
  p2 <- signif(v2_only/universe_size, 3)

  caption <- paste0("Proportion of DE significant ", v1_name, " genes : ", p1,
                    "\nProportion of DE significant ", v1v2_name, " genes : ", p1_2,
                    "\nProportion of DE significant ", v2_name, " genes : ", p2)

  plt <- invisible(plot(eulerr::euler(info4euler), fills = fills, quantities = quantities, ...))
  plt <- ggplotify::as.ggplot(plt)
  plt + ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
}

