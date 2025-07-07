#' Plot PCA Results for Gene Expression Data
#'
#' This function generates a PCA plot for a `SummarizedExperiment` object using precomputed PCA results.
#'
#' @param exp_data A `SummarizedExperiment` object containing gene expression data. PCA results must be stored in `exp_data@metadata` under the name specified by `res_name`.
#' @param color A character string specifying the column name in `colData(exp_data)` to use for coloring the points. Default is `NA`, meaning no coloring.
#' @param shape A character string specifying the column name in `colData(exp_data)` to use for shaping the points. Default is `NA`, meaning no shaping.
#' @param fname A character string specifying the file path to save the plot as a PDF. Default is `"results/qc/plot_PCA.pdf"`. Set to `NULL` to skip saving.
#' @param pcs An integer vector of length 2 specifying the principal components to plot on the x and y axes. Default is `c(1, 2)`.
#' @param res_name A character string specifying the name of the PCA results stored in `exp_data@metadata`. Default is `"pca_res"`.
#' @param id_name A character string specifying the column name in `colData(exp_data)` containing sample identifiers. Default is `"sample_id"`.
#' @param point_size A numeric value specifying the size of the points in the plot. Default is `2`.
#' @param ellipses if TRUE, `stat_ellipse` is added to the points, representing the 95% confidence-interval ellipse for each group in `color`.
#' @param out A character string indicating the output type: `"plotly"` (interactive Plotly plot) or `"ggplot"` (static ggplot). Default is `"plotly"`.
#'
#' @return A PCA plot, either as a `plotly` interactive object or a `ggplot` static object, depending on the `out` parameter.
#'
#' @details
#' The function creates a scatter plot of the specified principal components (`pcs[1]` and `pcs[2]`). Users can optionally:
#' - Color points based on a metadata column (`color`).
#' - Shape points based on a metadata column (`shape`).
#' - Save the static plot as a PDF if `fname` is not `NULL`.
#'
#' If PCA results are not precomputed using `pca_gexp`, the function will throw an error. The proportion of variance explained by the principal components is displayed on the axis labels.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_light guides guide_legend element_rect element_line element_text stat_ellipse
#' @importFrom plotly ggplotly
#' @importFrom dplyr left_join
#' @importFrom stringr str_c str_glue
#' @importFrom S4Vectors metadata
#' @importFrom rlang :=
#'
#' @export
#'
plot_pca <- function(exp_data,
                     color = NA,
                     shape = NA,
                     fname = "results/qc/plot_PCA.pdf",
                     pcs = c(1, 2),
                     res_name = "pca_res",
                     id_name = "sample_id",
                     ellipses = FALSE,
                     point_size = 2,
                     out = c("plotly", "ggplot")[1]) {

  if (is.null(metadata(exp_data)[[res_name]])) {
    stop("PCA results are missing. Please run `pca_gexp` before plotting.")
  }

  annotation <- as.data.frame(colData(exp_data))
  pca_res <- metadata(exp_data)[[res_name]]
  prop_var <- summary(pca_res)$importance[2,]
  prop_var_plt <- paste0(" (", format(prop_var * 100, digits = 2, trim = TRUE), "%)")

  plot_df <- pca_res$x[, pcs] %>% as.data.frame() %>% rownames_to_column(id_name)
  PCa <- paste0("PC", pcs[1])
  PCb <- paste0("PC", pcs[2])

  plot_df <- left_join(plot_df, annotation, by = id_name)

  # Check if color/shape exist in colData
  if (!is.na(color) && !color %in% colnames(annotation)) {
    warning(str_glue("Column '{color}' not found in colData. Disabling color aesthetic."))
    color <- NA
  }
  if (!is.na(shape) && !shape %in% colnames(annotation)) {
    warning(str_glue("Column '{shape}' not found in colData. Disabling shape aesthetic."))
    shape <- NA
  }
  pca_plot <- ggplot(plot_df, aes(!!sym(PCa), !!sym(PCb), label = !!sym(id_name)))

  if (is.na(color)) {
    pca_plot <- pca_plot +
      geom_point(aes(size = point_size))
  } else if (is.na(shape)) {
      pca_plot <- pca_plot +
        geom_point(aes(color = !!sym(color)), size = point_size)
  } else if (color != shape) {
    pca_plot <- pca_plot +
      geom_point(aes(color = !!sym(color),
                              shape = !!sym(shape)), size = point_size)
  } else {
    pca_plot <- pca_plot +
      geom_point(aes(color = !!sym(color),
                              shape = !!sym(color)), size = point_size)
  }

  pca_plot <- pca_plot +
    geom_hline(yintercept = 0, lty = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, lty = "dashed", color = "grey50") +
    labs(x = paste0("PC ", pcs[1], prop_var_plt[pcs[1]]),
                  y = paste0("PC ", pcs[2], prop_var_plt[pcs[2]])) +
    theme_light() +
    theme(panel.grid.minor = element_blank(),
                   panel.border = element_rect(color = "black"),
                   axis.ticks = element_line(color = "black"),
                   axis.text = element_text(color = "black"))

  # Legends
  # if (!is.na(color)) {
  #     pca_plot <- pca_plot + labs(color = color)
  # }
  if(ellipses) {
    if(!is.na(color)) {
      pca_plot <- pca_plot +
        stat_ellipse(aes(color = !!sym(color)))
    } else {
      pca_plot <- pca_plot + stat_ellipse()
    }
  }

  # Save plot if path provided
  if (!is.null(fname)) {
    ggsave(fname, pca_plot, device = "pdf")
  }

  if (out == "plotly") {
    return(ggplotly(pca_plot, tooltip = c(id_name, color)))
  } else {
    return(pca_plot)
  }
}
