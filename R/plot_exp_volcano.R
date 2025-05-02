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
  # Vérification des colonnes
  required_cols <- c("log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(diffexp))) {
    stop("Le data.frame diffexp doit contenir les colonnes : log2FoldChange, padj, gene")
  }

  # Calcul des catégories de régulation
  diffexp$Significance <- NA
  diffexp$Significance[diffexp$padj < 0.05 & diffexp$log2FoldChange > 1] <- "Upregulated"
  diffexp$Significance[diffexp$padj < 0.05 & diffexp$log2FoldChange < -1] <- "Downregulated"

  # Sélection des top 10 gènes (selon padj puis FC)
  top_genes <- dplyr::filter(diffexp, padj < 0.05 & abs(log2FoldChange) > 1)
  top_genes <- dplyr::arrange(top_genes, padj)
  top_genes <- dplyr::slice_head(top_genes, n = 10)

  # Ajouter les labels pour annotation
  diffexp$label <- ifelse(rownames(diffexp) %in% rownames(top_genes), rownames(diffexp), NA)

  # Définir la couleur des points : bleu pour upregulated, rouge pour downregulated
  diffexp$point_color <- "gray"  # Définir une couleur par défaut (gris pour les gènes non significatifs)
  diffexp$point_color[diffexp$Significance == "Upregulated"] <- "blue"
  diffexp$point_color[diffexp$Significance == "Downregulated"] <- "red"

  # Création du volcano plot
  plot <- ggplot2::ggplot(diffexp, ggplot2::aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = point_color  # Utiliser la couleur définie dans point_color
  )) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      color = "black",  # Texte des labels en noir
      max.overlaps = 100
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dotted",
      color = "darkgray"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1),
      linetype = "dotted",
      color = "darkgray"
    ) +
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

  return(plot)
}
