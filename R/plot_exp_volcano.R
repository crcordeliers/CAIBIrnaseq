#' Plot Volcano Plot of Differential Expression Results
#'
#' Generates a volcano plot based on the results of differential expression analysis,
#' highlighting upregulated and downregulated genes, with labels for top significant genes.
#'
#' @param diffexp A data frame containing differential expression results. Must include:
#'   - `log2FoldChange`: The log2 fold change values for each gene.
#'   - `padj`: The adjusted p-value for each gene.
#' @param nb The number of genes that have an annotation
#'
#' @return A `ggplot` object representing the volcano plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @importFrom forcats fct_rev
#' @export
#'
plot_exp_volcano <- function(diffexp, nb = 10, title = "Volcano Plot of Differential Expression") {
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

  # Sélection des top nb gènes
  top_genes <- diffexp |>
    dplyr::filter(Significance != "Not Significant") |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = nb)

  # Ajouter une colonne pour les labels
  diffexp$label <- ifelse(diffexp$gene %in% top_genes$gene, diffexp$gene, NA)

  # Palette de couleurs personnalisée
  color_palette <- c("Upregulated" = "#0072B2", "Downregulated" = "#D55E00", "Not Significant" = "gray80")

  # Définir les limites élargies pour dézoomer
  x_margin <- 2
  x_range <- range(diffexp$log2FoldChange, na.rm = TRUE)
  x_limits <- c(floor(x_range[1]) - x_margin, ceiling(x_range[2]) + x_margin)

  y_range <- -log10(diffexp$padj)
  y_limit <- ceiling(max(y_range, na.rm = TRUE)) + 1

  # Construction du plot
  vplot<- ggplot2::ggplot(diffexp, ggplot2::aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significance
  )) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 2.5,
      color = "black",
      max.overlaps = 100,
      box.padding = 0.4,
      max.iter = 10000,
      seed = 42
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "gray"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      color = "gray"
    ) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::scale_x_continuous(limits = x_limits) +
    ggplot2::scale_y_continuous(limits = c(0, y_limit)) +
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
  return (vplot)
}
