#' Plot Volcano Plot of Differential Expression Results
#'
#' Generates a volcano plot based on the results of differential expression analysis,
#' highlighting upregulated and downregulated genes, with labels for top significant genes.
#'
#' @param diffexp A data frame containing differential expression results. Must include:
#'   - `log2FoldChange`: The log2 fold change values for each gene.
#'   - `padj`: The adjusted p-value for each gene.
#' @param nb The number of genes that have an annotation
#' @param title An optional plot title. `NULL` (the default) omits the title.
#' @param lfc_threshold A single non-negative number giving the absolute log2 fold change
#' cutoff used to call a gene up-/down-regulated (i.e. `|log2FoldChange| > lfc_threshold`).
#' Default is `1` (a 2-fold change), a fixed, dataset-independent effect-size bar - not a
#' data-adaptive value like a percentile of the observed fold changes, since that would make
#' the same real effect size count as "significant" or not depending on the shape of each
#' dataset's fold-change distribution, breaking comparability across analyses.
#' @param padj_threshold A single number giving the adjusted p-value cutoff used to call a gene
#' significant. Default is `0.05`.
#' @param labelled_genes An optional character vector of gene names to label on the plot. 
#' If provided, these genes will be labeled instead of the top `nb` genes.
#' @param color_up Color used for upregulated genes (default is "#0072B2").
#' @param color_down Color used for downregulated genes (default is "#D55E00").
#' @param color_ns Color used for non-significant genes (default is "gray80").
#' @param fname Optional file name (including path) to save the static ggplot version of the volcano plot.
#' @param out A character string indicating the output type: `"plotly"` (interactive Plotly plot)
#' or `"ggplot"` (static ggplot). Default is `"plotly"`. Only affects the returned object; when
#' `fname` is provided, the static ggplot version is always what gets saved to file.
#'
#' @return A volcano plot, either as a `ggplot` static object or a `plotly` interactive object,
#' depending on the `out` parameter.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_hline geom_vline labs theme_minimal element_text ggsave
#' @importFrom dplyr filter arrange slice_head
#' @importFrom forcats fct_rev
#' @importFrom fs path_dir
#' @importFrom plotly ggplotly
#' @importFrom ggrepel geom_text_repel
#' @export
#'
plot_exp_volcano <- function(diffexp, nb = 10,
                             title = NULL,
                             lfc_threshold = 1,
                             padj_threshold = 0.05,
                             labelled_genes = NULL,
                             color_up = "#0072B2",
                             color_down = "#D32F2F",
                             color_ns = "gray80",
                             fname = NULL,
                             out = c("plotly", "ggplot")[1]) {
  # Check required columns
  required_cols <- c("log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(diffexp))) {
    stop("diffexp must contain columns: log2FoldChange and padj.")
  }
  if (length(lfc_threshold) != 1 || lfc_threshold < 0) {
    stop("`lfc_threshold` must be a single non-negative number.")
  }

  # Handle gene names
  if (!"gene" %in% colnames(diffexp)) {
    diffexp$gene <- rownames(diffexp)
  }

  lfc_low <- -lfc_threshold
  lfc_high <- lfc_threshold

  # Define significance groups
  diffexp$Significance <- case_when(
    diffexp$padj < padj_threshold & diffexp$log2FoldChange > lfc_high ~ "Upregulated",
    diffexp$padj < padj_threshold & diffexp$log2FoldChange < lfc_low ~ "Downregulated",
    TRUE ~ "Not Significant"
  )

  if(is.null(labelled_genes)) {
    # Select top nb genes
    top_genes <- diffexp |>
      filter(Significance != "Not Significant") |>
      arrange(padj) |>
      slice_head(n = nb)

    top_genes <-  top_genes$gene
  } else {
    # Use provided labelled_genes
    top_genes <- intersect(labelled_genes, diffexp$gene)
  }

  # Label shown on-plot (via geom_text)
  diffexp$repel_label <- ifelse(diffexp$gene %in% top_genes, diffexp$gene, NA)

  # Custom plotly hover text: gene name, then log2FC and padj rounded to 3 significant figures.
  diffexp$hover_text <- paste0(
    "Gene: ", diffexp$gene,
    "<br>log2FC: ", signif(diffexp$log2FoldChange, 3),
    "<br>padj: ", signif(diffexp$padj, 3)
  )

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
    color = Significance,
    text = hover_text
  )) +
    geom_point(alpha = 0.7, size = 1.5)

  if (out == "ggplot") {
    vplot <- vplot +
      geom_text_repel(
        aes(label = repel_label),
        size = 2.5,
        color = "black",
        max.overlaps = Inf
      )
  } else {
    vplot <- vplot +
      geom_text(
        aes(label = repel_label),
        size = 2.5,
        color = "black",
        vjust = -0.8,
        nudge_y = y_limit * 0.03
      )
  }

  vplot <- vplot +
    geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = "dashed",
      color = "gray"
    ) +
    geom_vline(
      xintercept = c(lfc_low, lfc_high),
      linetype = "dashed",
      color = "gray"
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = c(0, y_limit)) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 Adjusted p-value",
      color = "Regulation"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  if(!is.null(fname)) {
    dir.create(path_dir(fname), recursive = TRUE, showWarnings = FALSE)
    ggsave(fname, vplot)
  }

  if (out == "plotly") {
    return(ggplotly(vplot, tooltip = "text"))
  } else {
    return(vplot)
  }
}
