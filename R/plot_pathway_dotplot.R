#' Plot Pathway Analysis Results as a Dot Plot
#'
#' Creates a dot plot to visualize pathway enrichment analysis results, supporting both ORA and FGSEA formats.
#'
#' @param exp_data A `SummarizedExperiment` object containing pathway analysis results stored in the `metadata()`.
#' @param score_name A character string indicating the metadata field where pathway results are stored. Default is `"resultsORA"`.
#' @param top_n Number of top pathways to display. Default is `10`.
#' @param maxPval Maximum adjusted p-value for filtering pathways (if `usePval = TRUE`). Default is `0.05`.
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
                                 score_name,
                                 top_n = 10,
                                 maxPval = 0.05) {

  results <- metadata(exp_data)[[score_name]]

  results_type <- case_when(
    "geneHits" %in% colnames(results) ~ "ORA",
    "leadingEdge" %in% colnames(results) ~ "GSEA",
    TRUE ~ NA_character_
  )

  if (results_type == "ORA") {
    message("Plotting ORA-type results...")

    results <- results |>
      arrange(padj) |>
      filter(padj < maxPval) |>
      slice(1:top_n) |>
      rowwise() |>
      mutate(geneRatioNum = .fraction_to_decimal(geneRatio)) |>
      ungroup()

    results <- results %>%
      mutate(
        pathway = factor(pathway, levels = rev(pathway)),
        logpadj = -log10(padj)
      )

    p <- ggplot(results, aes(x = logpadj, y = pathway, size = geneRatioNum, fill = logpadj)) +
      scale_x_continuous(
        name = "Adjusted p-value",
        labels = function(x) signif(10^-x, 2),
        breaks = scales::pretty_breaks(n = 5)
      )

  } else if (results_type == "GSEA") {
    message("Plotting GSEA-type results...")

    results <- results %>%
      arrange(padj) |>
      filter(padj < maxPval) %>%
      slice(1:top_n) %>%
      arrange(NES) %>%
      mutate(
        geneRatioNum = lengths(leadingEdge) / size,
        pathway = factor(pathway, levels = pathway),
        logpadj = -log10(padj))

    p <- ggplot(results, aes(x = NES, y = pathway, size = geneRatioNum, fill = logpadj)) +
      scale_x_continuous(
        name = "Normalized Enrichment Score (NES)",
        expand = c(0.1, 0.1)
      ) +
      geom_vline(xintercept = 0, linetype = "dashed")

  } else {
    stop("Unsupported results type format: expected columns `geneRatio` or `leadingEdge`.")
  }

  p <- p +
    scale_fill_distiller(
      name = "Adjusted p-value",
      palette = "Reds", direction = 1,
      labels = function(x) signif(10^-x, 2),
      limits = c(-log10(maxPval), max(results$logpadj, na.rm = TRUE))
    ) +
    geom_segment(aes(xend = 0, yend = pathway), size = 0.5, color = "grey", linetype = "dotted") +
    geom_point(pch = 21, color = "black") +
    scale_size_continuous(range = c(2, 10)) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "), width = 30)) +
    ggtitle("Pathway Analysis") +
    labs(y = "Pathway", size = "Gene Ratio") +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = "gray90")
    )

  return(p)
}

#' @importFrom stringr str_split
.fraction_to_decimal <- function(fraction) {
  num_denum <- str_split_fixed(fraction, "/", 2) |>
    as.numeric()
  return(num_denum[1]/num_denum[2])
}

