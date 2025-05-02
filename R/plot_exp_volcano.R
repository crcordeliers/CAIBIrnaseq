plot_exp_volcano <- function(diffexp, title = "Volcano Plot of Differential Expression") {
  # Vérification des colonnes requises
  required_cols <- c("log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(diffexp))) {
    stop("Le data.frame diffexp doit contenir les colonnes : log2FoldChange et padj.")
  }

  # Gestion du nom des gènes
  if (!"gene" %in% colnames(diffexp)) {
    diffexp$gene <- rownames(diffexp)
  }

  # Définir les groupes de significativité
  diffexp$Significance <- dplyr::case_when(
    diffexp$padj < 0.05 & diffexp$log2FoldChange > 1 ~ "Upregulated",
    diffexp$padj < 0.05 & diffexp$log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )

  # Sélection des top 10 gènes
  top_genes <- diffexp |>
    dplyr::filter(Significance != "Not Significant") |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = 10)

  # Ajouter une colonne pour les labels
  diffexp$label <- ifelse(diffexp$gene %in% top_genes$gene, diffexp$gene, NA)

  # Palette de couleurs personnalisée
  color_palette <- c("Upregulated" = "#0072B2", "Downregulated" = "#D55E00", "Not Significant" = "gray80")

  # Construction du plot
  ggplot2::ggplot(diffexp, ggplot2::aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significance
  )) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 3,
      max.overlaps = 100,
      box.padding = 0.4,
      max.iter = 10000,
      seed = 42
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "gray40"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      color = "gray40"
    ) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted p-value",
      color = "Regulation"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}
