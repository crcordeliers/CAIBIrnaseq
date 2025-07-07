#' Plot Volcano Plot of Differential Expression Results
#'
#' Generates a volcano plot based on the results of differential expression analysis,
#' highlighting upregulated and downregulated genes, with labels for top significant genes.
#'
#' @param diffexp A data frame containing differential expression results. Must include:
#'   - `log2FoldChange`: The log2 fold change values for each gene.
#'   - `padj`: The adjusted p-value for each gene.
#' @param nb The number of genes that have an annotation
#' @param color_up Color used for upregulated genes (default is "#0072B2").
#' @param color_down Color used for downregulated genes (default is "#D55E00").
#' @param color_ns Color used for non-significant genes (default is "gray80").
#'
#' @return A `ggplot` object representing the volcano plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @importFrom forcats fct_rev
#' @export
#'
plot_exp_volcano <- function(diffexp, nb = 10, title = "Volcano Plot of Differential Expression",
                             color_up = "#0072B2",
                             color_down = "#D55E00",
                             color_ns = "gray80") {
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
  diffexp$Significance <- case_when(
    diffexp$padj < 0.05 & diffexp$log2FoldChange > 1 ~ "Upregulated",
    diffexp$padj < 0.05 & diffexp$log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )

  # Sélection des top nb gènes
  top_genes <- diffexp |>
    filter(Significance != "Not Significant") |>
    arrange(padj) |>
    slice_head(n = nb)

  # Ajouter une colonne pour les labels
  diffexp$label <- ifelse(diffexp$gene %in% top_genes$gene, diffexp$gene, NA)

  # Palette personnalisée depuis les paramètres
  color_palette <- c(
    "Upregulated" = color_up,
    "Downregulated" = color_down,
    "Not Significant" = color_ns
  )

  # Définir les limites élargies
  x_margin <- 2
  x_range <- range(diffexp$log2FoldChange, na.rm = TRUE)
  x_limits <- c(floor(x_range[1]) - x_margin, ceiling(x_range[2]) + x_margin)

  y_range <- -log10(diffexp$padj)
  y_limit <- ceiling(max(y_range, na.rm = TRUE)) + 1

  # Construction du graphique
  vplot <- ggplot(diffexp, aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significance
  )) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_text_repel(
      aes(label = label),
      size = 2.5,
      color = "black",
      max.overlaps = 100,
      box.padding = 0.4,
      max.iter = 10000,
      seed = 42
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "gray"
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      color = "gray"
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = c(0, y_limit)) +
    labs(
      title = title,
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~adjusted~italic(p)-value),
      color = "Regulation"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  return(vplot)
}
