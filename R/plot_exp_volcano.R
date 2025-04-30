#' Plot Volcano Plot of Differential Expression Results
#'
#' Generates an interactive volcano plot based on the results of differential expression analysis, highlighting upregulated and downregulated genes.
#'
#' @param diffexp A data frame containing differential expression results. Must include the following columns:
#'   - `log2FoldChange`: The log2 fold change values for each gene.
#'   - `padj`: The adjusted p-value for each gene.
#'
#' @return A `plotly` object representing the interactive volcano plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom plotly ggplotly
#' @importFrom dplyr filter arrange slice_head
#' @importFrom forcats fct_rev
#' @export
#'
plot_exp_volcano <- function(diffexp) {
  ## --------------------------- ##
  ## Verifications / Validation  ##
  ## --------------------------- ##
  if (!is.data.frame(diffexp)) {
    stop("diffexp must be a data.frame.")
  }

  required_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(diffexp))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from diffexp: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(diffexp$log2FoldChange)) {
    stop("'log2FoldChange' column must be numeric.")
  }
  if (!is.numeric(diffexp$padj)) {
    stop("'padj' column must be numeric.")
  }

  if (is.null(rownames(diffexp))) {
    stop("diffexp must have rownames corresponding to gene names.")
  }

  ## --------------------------- ##
  ##         Main Function        ##
  ## --------------------------- ##

  diffexp$Significance <- NA
  diffexp$Significance[diffexp$padj < 0.05 & diffexp$log2FoldChange < -1] <- "Downregulated"
  diffexp$Significance[diffexp$padj < 0.05 & diffexp$log2FoldChange > 1] <- "Upregulated"

  top_genes <- diffexp |>
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = 20)

  volcanoPlot <- ggplot2::ggplot(
    diffexp,
    ggplot2::aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = forcats::fct_rev(Significance),
      text = paste(
        "Gene:", rownames(diffexp),
        "<br>Log2FC:", signif(log2FoldChange, 3),
        "<br>-log10(padj):", signif(-log10(padj), 3)
      )
    )
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "darkgray") +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "darkgray") +
    ggplot2::labs(
      title = "Volcano Plot of Differential Expression",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted p-value",
      color = "Significance"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2:plot_exp_volcano <- function(diffexp) {
      ## --------------------------- ##
      ## Verifications / Validation  ##
      ## --------------------------- ##
      if (!is.data.frame(diffexp)) {
        stop("diffexp must be a data.frame.")
      }

      required_cols <- c("log2FoldChange", "padj")
      missing_cols <- setdiff(required_cols, colnames(diffexp))
      if (length(missing_cols) > 0) {
        stop("The following required columns are missing from diffexp: ", paste(missing_cols, collapse = ", "))
      }

      if (!is.numeric(diffexp$log2FoldChange)) {
        stop("'log2FoldChange' column must be numeric.")
      }
      if (!is.numeric(diffexp$padj)) {
        stop("'padj' column must be numeric.")
      }

      if (is.null(rownames(diffexp))) {
        stop("diffexp must have rownames corresponding to gene names.")
      }

      ## --------------------------- ##
      ##         Main Function        ##
      ## --------------------------- ##

      diffexp$Significance <- NA
      diffexp$Significance[diffexp$padj < 0.05 & diffexp$log2FoldChange < -1] <- "Downregulated"
      diffexp$Significance[diffexp$padj < 0.05 & diffexp$log2FoldChange > 1] <- "Upregulated"

      top_genes <- diffexp |>
        dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) |>
        dplyr::arrange(padj) |>
        dplyr::slice_head(n = 20)

      volcanoPlot <- ggplot2::ggplot(
        diffexp,
        ggplot2::aes(
          x = log2FoldChange,
          y = -log10(padj),
          color = forcats::fct_rev(Significance),
          text = paste(
            "Gene:", rownames(diffexp),
            "<br>Log2FC:", signif(log2FoldChange, 3),
            "<br>-log10(padj):", signif(-log10(padj), 3)
          )
        )
      ) +
        ggplot2::geom_point(alpha = 0.6, size = 2) +
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

      # interactivePlot <- plotly::ggplotly(volcanoPlot, tooltip = "text")

      return(volcanoPlot)
    }
  theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  # interactivePlot <- plotly::ggplotly(volcanoPlot, tooltip = "text")

  return(volcanoPlot)
}
