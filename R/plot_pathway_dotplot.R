#' Plot Pathway Analysis Results as a Dot Plot
#'
#' Creates a dot plot to visualize pathway analysis results, supporting both ORA and FGSEA formats.
#'
#' @param exp_data A `SummarizedExperiment` object containing experimental data with pathway analysis results in `metadata()`.
#' @param score_name A character string indicating the metadata field where results are stored. Default is `"resultsORA"`.
#' @param top_n An integer for the number of top pathways to plot. Default is `10`.
#' @param maxPval A numeric for the maximum adjusted p-value to include pathways. Default is `0.05`.
#'
#' @return A `ggplot2` dot plot showing pathway enrichment results.
#'
#' @details
#' The function supports two result formats:
#' - **ORA** (Over-Representation Analysis): expects columns `PAdj`, `GeneRatio`, and `Pathway`.
#' - **FGSEA**: expects columns `padj`, `size`, and `pathway`.
#'
#' The x-axis shows the -log10 of adjusted p-values, y-axis lists pathways, dot size shows gene ratio.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous scale_color_gradient scale_size_continuous scale_y_discrete ggtitle labs theme element_text element_blank element_line
#' @importFrom dplyr arrange filter slice_head mutate
#' @importFrom scales pretty_breaks
#' @importFrom stringr str_wrap str_replace_all
#' @importFrom S4Vectors metadata
#'
#' @export
#'
plot_pathway_dotplot <- function(exp_data, score_name = "resultsORA", top_n = 10, maxPval = 0.05, usePval = TRUE) {
  results <- S4Vectors::metadata(exp_data)[[score_name]]

  if ("PValue" %in% colnames(results)) {
    # ORA format
    if (usePval) {
      results <- results |>
        dplyr::arrange(PAdj) |>
        dplyr::filter(PAdj < maxPval) |>
        dplyr::slice_head(n = top_n)
    } else {
      # Pas de filtrage ni tri par p-value, on trie par GeneRatio décroissant
      results <- results |>
        dplyr::mutate(
          GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))
        ) |>
        dplyr::arrange(desc(GeneRatioNum)) |>
        dplyr::slice_head(n = top_n)
    }
    results <- results |>
      dplyr::mutate(
        GeneRatioNum = ifelse(exists("GeneRatioNum"), GeneRatioNum,
                              as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))),
        Pathway = factor(Pathway, levels = rev(Pathway)),
        Size = GeneRatioNum,
        logpadj = ifelse(usePval, -log10(PAdj), NA_real_)
      )
  } else if ("pval" %in% colnames(results)) {
    # FGSEA format
    if (usePval) {
      results <- results |>
        dplyr::arrange(padj) |>
        dplyr::filter(padj < maxPval) |>
        dplyr::slice_head(n = top_n)
    } else {
      # Pas de filtrage ni tri par p-value, on trie par size décroissant
      results <- results |>
        dplyr::arrange(desc(size)) |>
        dplyr::slice_head(n = top_n)
    }
    results <- results |>
      dplyr::mutate(
        GeneRatioNum = ifelse(usePval, size / max(size), size / max(size)),
        Pathway = factor(pathway, levels = rev(pathway)),
        Size = GeneRatioNum,
        logpadj = ifelse(usePval, -log10(padj), NA_real_)
      )
  } else {
    stop("Unsupported results format: expected columns `PAdj` or `padj`.")
  }

  if (usePval) {
    p <- ggplot2::ggplot(results, ggplot2::aes(x = logpadj, y = Pathway, size = Size, color = logpadj)) +
      ggplot2::scale_x_continuous(labels = function(x) 10^-x, breaks = scales::pretty_breaks(n = 5)) +
      ggplot2::scale_color_gradient(
        low = "#df6664", high = "#387eb9",
        labels = function(x) 10^-x,
        breaks = scales::pretty_breaks(n = 6),
        limits = c(-log10(maxPval), max(results$logpadj, na.rm = TRUE))
      ) +
      ggplot2::labs(x = "Adjusted p-value", color = "Adjusted p-value")
  } else {
    p <- ggplot2::ggplot(results, ggplot2::aes(x = Size, y = Pathway, size = Size)) +
      ggplot2::geom_point(color = "#387eb9") +
      ggplot2::labs(x = "Gene Ratio", color = NULL)
  }

  p +
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
}

