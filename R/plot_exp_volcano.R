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
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal element_text ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @importFrom forcats fct_rev
#' @importFrom fs path_dir
#' @export
#'
plot_exp_volcano <- function(diffexp, nb = 10, title = "Volcano Plot of Differential Expression",
                             color_up = "#a8250eff",
                             color_down = "#2954b1ff",
                             color_ns = "gray80",
                             fname = NULL) {
  # Check required columns
  required_cols <- c("log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(diffexp))) {
    stop("diffexp must contain columns: log2FoldChange and padj.")
  }

  # Handle gene names
  if (!"gene" %in% colnames(diffexp)) {
    diffexp$gene <- rownames(diffexp)
  }

  # Define significance groups
  diffexp$Significance <- case_when(
    diffexp$padj < 0.05 & diffexp$log2FoldChange > 1 ~ "Upregulated",
    diffexp$padj < 0.05 & diffexp$log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )

  # Select top nb genes
  top_genes <- diffexp |>
    filter(Significance != "Not Significant") |>
    arrange(padj) |>
    slice_head(n = nb)

  # Add label column
  diffexp$label <- ifelse(diffexp$gene %in% top_genes$gene, diffexp$gene, NA)

  # Custom color palette
  color_palette <- c(
    "Upregulated" = color_up,
    "Downregulated" = color_down,
    "Not Significant" = color_ns
  )

  # Define expanded axis limits
  x_margin <- 2
  x_range <- range(diffexp$log2FoldChange, na.rm = TRUE)
  x_limits <- c(floor(x_range[1]) - x_margin, ceiling(x_range[2]) + x_margin)

  y_range <- -log10(diffexp$padj)
  y_limit <- ceiling(max(y_range, na.rm = TRUE)) + 1

  # Build plot
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

  if(!is.null(fname)) {
    dir.create(path_dir(fname), recursive = TRUE, showWarnings = FALSE)
    ggsave(fname, vplot)
  }

  return(vplot)
}
