#' PCA Plot for Gene Expression Data
#'
#' @param exp_data A SummarizedExperiment object containing gene expression data.
#' @param color Variable in the data for coloring the points (default is NA).
#' @param shape Variable in the data for changing the shape of the points (default is NA).
#' @param fname File name to save the plot (default is "results/qc/plot_PCA.pdf").
#' @param pcs A vector of two integers specifying which principal components to plot (default is c(1, 2)).
#' @param res_name Name of the PCA result stored in the metadata of the `exp_data` object (default is "pca_res").
#' @param id_name Column name for sample identifiers (default is "sample_id").
#' @param point_size Size of the points in the plot (default is 2).
#' @param color_in Either "color" or "fill" specifying how to handle color mapping (default is "color").
#' @param out Output format: either "plotly" for interactive or "ggplot" for static (default is "plotly").
#'
#' @returns A ggplot or plotly object depending on the output format selected.
#' @export
#' @importFrom ggplot2 theme element_blank element_rect element_line element_text guides guide_legend
#'
plot_pca <- function(exp_data,
                     color = NA,
                     shape = NA,
                     fname = "results/qc/plot_PCA.pdf",
                     pcs = c(1, 2),
                     res_name = "pca_res",
                     id_name = "sample_id",
                     point_size = 2,
                     color_in = c("color", "fill")[1],
                     out = c("plotly", "ggplot")[1]) {

  # Ensure PCA results exist
  if (is.null(exp_data@metadata[[res_name]])) {
    stop("PCA results are missing. Please run `pca_gexp` before plotting.")
  }

  # Get PCA results and sample annotations
  annotation <- as.data.frame(colData(exp_data))
  pca_res <- exp_data@metadata$pca_res

  # Calculate proportion of variance explained by each principal component
  prop_var <- summary(pca_res)$importance[2,]
  prop_var_plt <- paste0(" (", format(prop_var * 100, digits = 2, trim = TRUE), "%)")

  # Prepare data for plotting (PCA components and annotations)
  plot_df <- pca_res$x[, pcs] %>% as.data.frame() %>% rownames_to_column(id_name)
  PCa <- paste0("PC", pcs[1])
  PCb <- paste0("PC", pcs[2])

  # Join PCA data with sample annotations
  plot_df <- left_join(plot_df, annotation, by = id_name)

  # Create ggplot object
  pca_plot <- ggplot(plot_df, aes(label = !!sym(id_name), x = !!sym(PCa), y = !!sym(PCb))) +
    geom_point(size = point_size)

  # Add color and shape aesthetics
  if (!is.na(color)) {
    pca_plot <- pca_plot + aes(color = !!sym(color))
  }
  if (!is.na(shape)) {
    pca_plot <- pca_plot + aes(shape = !!sym(shape))
  }

  # Style adjustments for the plot
  pca_plot <- pca_plot +
    geom_hline(yintercept = 0, lty = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, lty = "dashed", color = "grey50") +
    labs(x = paste0("PC ", pcs[1], prop_var_plt[pcs[1]]),
         y = paste0("PC ", pcs[2], prop_var_plt[pcs[2]])) +
    theme_light() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black"),
          axis.ticks = ggplot2::element_line(color = "black"),
          axis.text = ggplot2::element_text(color = "black"))

  # Handle color or fill legends
  if (color_in == "color") {
    pca_plot <- pca_plot + labs(color = color)
  } else {
    pca_plot <- pca_plot + labs(fill = color) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21)))
  }

  # Convert ggplot to plotly for interactive plotting
  plt_pca_int <- ggplotly(pca_plot, tooltip = c(id_name, color))

  # Save the plot if a file name is provided
  if (!is.null(fname)) {
    ggsave(fname, pca_plot, create.dir = TRUE)
  }

  # Return either a static ggplot or interactive plotly object
  if (out == "plotly") {
    return(plt_pca_int)
  } else {
    return(pca_plot)
  }
}
