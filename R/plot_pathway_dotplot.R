#' Plot Pathway Analysis Results as a Dot Plot
#'
#' Creates a dot plot to visualize pathway enrichment analysis results, supporting both ORA and FGSEA formats.
#'
#' @param exp_data A `SummarizedExperiment` object containing pathway analysis results stored in the `metadata()`.
#' @param score_name A character string indicating the metadata field where pathway results are stored. Default is `"resultsORA"`.
#' @param top_n Number of top pathways to display. Default is `10`.
#' @param maxPval Maximum adjusted p-value for filtering pathways (if `usePval = TRUE`). Default is `0.05`.
#' @param usePval Logical; whether to filter and rank by adjusted p-values (`TRUE`) or by gene ratio / size (`FALSE`). Default is `TRUE`.
#'
#' @return A `ggplot2` dot plot object showing pathway enrichment.
#'
#' @details
#' The function supports two types of enrichment results:
#' \itemize{
#'   \item{ORA (Over-Representation Analysis): requires columns `PAdj`, `GeneRatio`, and `Pathway`}
#'   \item{FGSEA: requires columns `padj`, `size`, and `pathway`}
#' }
#'
#' When `usePval = TRUE`, the plot will show -log10 adjusted p-values on the x-axis, colored by significance.
#' When `usePval = FALSE`, the plot will rank and size pathways by gene ratio or enrichment size.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous scale_color_gradient scale_size_continuous scale_y_discrete ggtitle labs theme element_text element_blank element_line
#' @importFrom dplyr arrange filter slice_head mutate
#' @importFrom scales pretty_breaks
#' @importFrom stringr str_wrap str_replace_all
#' @importFrom S4Vectors metadata
#'
#' @export
#'
plot_pathway_dotplot <- function(exp_data,
                                 score_name = "resultsORA",
                                 top_n = 10,
                                 maxPval = 0.05,
                                 usePval = TRUE) {

  results <- S4Vectors::metadata(exp_data)[[score_name]]

  if ("PValue" %in% colnames(results)) {
    # ORA format
    if (usePval) {
      results <- results |>
        dplyr::arrange(PAdj) |>
        dplyr::filter(PAdj < maxPval) |>
        dplyr::slice_head(n = top_n)
    } else {
      results <- results |>
        dplyr::mutate(GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))) |>
        dplyr::arrange(desc(GeneRatioNum)) |>
        dplyr::slice_head(n = top_n)
    }

    results <- results |>
      dplyr::mutate(
        GeneRatioNum = if (!"GeneRatioNum" %in% colnames(.)) as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)) else GeneRatioNum,
        Pathway = factor(Pathway, levels = rev(Pathway)),
        Size = GeneRatioNum,
        logpadj = if (usePval) -log10(PAdj) else NA_real_
      )

  } else if ("pval" %in% colnames(results)) {
    # FGSEA format
    if (usePval) {
      results <- results |>
        dplyr::arrange(padj) |>
        dplyr::filter(padj < maxPval) |>
        dplyr::slice_head(n = top_n)
    } else {
      results <- results |>
        dplyr::arrange(desc(size)) |>
        dplyr::slice_head(n = top_n)
    }

    results <- results |>
      dplyr::mutate(
        GeneRatioNum = size / max(size),
        Pathway = factor(pathway, levels = rev(pathway)),
        Size = GeneRatioNum,
        logpadj = if (usePval) -log10(padj) else NA_real_
      )

  } else {
    stop("Unsupported results format: expected columns `PAdj` or `padj`.")
  }

  if (usePval) {
    p <- ggplot2::ggplot(results, ggplot2::aes(x = logpadj, y = Pathway, size = Size, color = logpadj)) +
      ggplot2::scale_x_continuous(
        name = "Adjusted p-value",
        labels = function(x) signif(10^-x, 2),
        breaks = scales::pretty_breaks(n = 5)
      ) +
      ggplot2::scale_color_gradient(
        name = "Adjusted p-value",
        low = "#df6664", high = "#387eb9",
        labels = function(x) signif(10^-x, 2),
        breaks = scales::pretty_breaks(n = 6),
        limits = c(-log10(maxPval), max(results$logpadj, na.rm = TRUE))
      )
  } else {
    p <- ggplot2::ggplot(results, ggplot2::aes(x = Size, y = Pathway, size = Size)) +
      ggplot2::geom_point(color = "#387eb9") +
      ggplot2::labs(x = "Gene Ratio", color = NULL)
  }

  p <- p +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(range = c(2, 10)) +
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "), width = 30)) +
    ggplot2::ggtitle("Pathway Analysis") +
    ggplot2::labs(y = "Pathway", size = "Gene Ratio") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(size = 10, angle = ifelse(usePval, 45, 0), hjust = 1),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "gray90")
    )

  return(p)
}
