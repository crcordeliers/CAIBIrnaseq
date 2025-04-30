#' Plot Volcano Plot of Differential Expression Results
#'
#' Generates a volcano plot based on the results of differential expression analysis,
#' highlighting upregulated and downregulated genes, with labels for top significant genes.
#'
#' @param diffexp A data frame containing differential expression results. Must include:
#'   - `log2FoldChange`: The log2 fold change values for each gene.
#'   - `padj`: The adjusted p-value for each gene.
#'
#' @return A `ggplot` object representing the volcano plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @importFrom forcats fct_rev
#' @export
#'
plot_exp_volcano <- function(diffexp) {
  # Ajouter une colonne de signification
  diffexp$Significance <- NA
  diffexp$Significance[diffexp$padj < 0.05 & diffexp$log2FoldChange < -1] <- "Downregulated"
  diffexp$Significance[diffexp$padj < 0.05 & diffexp$log2FoldChange > 1] <- "Upregulated"

  # Sélection des top 10 gènes significatifs
  top_genes <- dplyr::filter(diffexp, padj < 0.05 & abs(log2FoldChange) > 1) |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = 10)

  # Ajoute une colonne 'label' pour les gènes à annoter
  diffexp$label <- ifelse(diffexp$gene %in% top_genes$gene, diffexp$gene, NA)

  # Création du volcano plot
  volcanoPlot <- ggplot2::ggplot(diffexp, ggplot2::aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = forcats::fct_rev(Significance)
  )) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), max.overlaps = 100) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "darkgray") +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "darkgray") +
    ggplot2::labs(
      title = "Volcano Plot of Differential Expression",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted p-value",
      color = "Significance"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  return(volcanoPlot)
}
