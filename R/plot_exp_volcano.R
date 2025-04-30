#' Volcano Plot for DESeq2 Results
#'
#' Creates a volcano plot from DESeq2 differential expression results.
#'
#' @param resLFC A data.frame or tibble containing the DESeq2 results with `log2FoldChange`, `padj`, and gene annotations.
#' @param contrast_name Character. Name of the contrast shown on the x-axis label.
#' @param lfc_tresh Numeric. Log2 fold-change threshold to define up/down-regulated genes. Default is 2.5.
#' @param padj_tresh Numeric. Adjusted p-value threshold for significance. Default is 0.05.
#' @param ngenes_show Integer. Number of top significant genes to label (split equally up/down). Default is 10.
#' @param gene_col Character. Column name in `resLFC` containing gene names. Default is "external_gene_name".
#' @param annot_pos Numeric vector of length 2. Position to place DEG annotation text (x, y). Use `Inf` or `-Inf` for corners.
#'
#' @return A `ggplot2` object displaying the volcano plot.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual annotate labs guides theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter slice_min pull case_when mutate
#' @importFrom rlang .data
#' @importFrom rlang enquo quo_name
#'
plot_exp_volcano <- function(resLFC,
                              contrast_name = "",
                              lfc_tresh = 2.5,
                              padj_tresh = 0.05,
                              ngenes_show = 10,
                              gene_col = "external_gene_name",
                              annot_pos = c(Inf, Inf)) {

  gene_sym <- rlang::ensym(gene_col)

  # Select top upregulated and downregulated genes
  genes_up <- resLFC |>
    dplyr::filter(.data$log2FoldChange > lfc_tresh) |>
    dplyr::slice_min(order_by = .data$padj, n = ngenes_show/2) |>
    dplyr::pull(!!gene_sym)

  genes_down <- resLFC |>
    dplyr::filter(.data$log2FoldChange < -lfc_tresh) |>
    dplyr::slice_min(order_by = .data$padj, n = ngenes_show/2) |>
    dplyr::pull(!!gene_sym)

  genes_label <- c(genes_up, genes_down)

  fc_leg <- paste("FC", contrast_name)
  ndegs <- sum(resLFC$padj < padj_tresh, na.rm = TRUE)
  deg_annot <- paste0("atop(bold('DEGs:'), ", ndegs, "/", nrow(resLFC), ")")

  # Prepare data for plot
  data2plot <- resLFC |>
    dplyr::mutate(
      lpval = -log10(.data$padj),
      signif = dplyr::case_when(
        .data$log2FoldChange > lfc_tresh & .data$padj < padj_tresh ~ "up",
        .data$log2FoldChange < -lfc_tresh & .data$padj < padj_tresh ~ "down",
        .data$padj < padj_tresh ~ "signif",
        TRUE ~ "ns"
      ),
      gene_name = .data[[gene_col]],
      label = ifelse(gene_name %in% genes_label, gene_name, NA)
    )

  # Define annotation positioning
  hjust_value <- dplyr::case_when(
    annot_pos[1] == Inf ~ 1,
    annot_pos[1] == -Inf ~ -1,
    TRUE ~ 0
  )

  vjust_value <- dplyr::case_when(
    annot_pos[2] == Inf ~ 1,
    annot_pos[2] == -Inf ~ -1,
    TRUE ~ 0
  )

  # Create plot
  plt <- ggplot2::ggplot(data2plot, ggplot2::aes(.data$log2FoldChange, lpval)) +
    ggplot2::geom_point(ggplot2::aes(color = signif), alpha = 0.7) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3) +
    ggplot2::annotate("text", x = annot_pos[1], y = annot_pos[2],
                      vjust = vjust_value, hjust = hjust_value,
                      label = deg_annot, parse = TRUE, size = 3.5) +
    ggplot2::scale_color_manual(values = c("up" = "#ef8a62",
                                           "down" = "#67a9cf",
                                           "signif" = "grey60",
                                           "ns" = "grey90")) +
    ggplot2::labs(x = bquote("log"[2]*.(fc_leg)), y = expression(-log[10]*" adj. p-value")) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_minimal()

  # Return the plot object
  return(plt)
}
