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
  volcanoPlot <- ggplot2::ggplot(data2plot, ggplot2::aes(log2FoldChange, lpval)) +
    ggplot2::geom_point(ggplot2::aes(color = signif), alpha = 0.7) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3) +
    ggplot2::annotate("text", x = annot_pos[1], y = annot_pos[2],
                      vjust = vjust_value, hjust = hjust_value,
                      label = deg_annot, size = 3.5, parse = TRUE) +
    ggplot2::scale_color_manual(values = c("up" = "#ef8a62",
                                           "down" = "#67a9cf",
                                           "signif" = "grey60",
                                           "ns" = "grey90")) +
    ggplot2::labs(x = bquote("log"[2]*.(fc_leg)), y = expression(-"log"[10]*"adj. p-value")) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_minimal()

  return(volcanoPlot)
}
