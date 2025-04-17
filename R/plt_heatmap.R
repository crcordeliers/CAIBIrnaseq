#' Plot Heatmap with Optional Annotations and Customization
#'
#' This function generates a heatmap from a data frame or matrix of values, with the option to include
#' additional annotations and customize various aspects of the plot such as fonts, track colors, and size.
#'
#' @param data4hm A list containing the data for the heatmap with the several elements
#' @param annotations A data frame or list of annotations to overlay on the heatmap. Default is `NA`.
#' @param fontsize The font size for text in the heatmap. Default is 9.
#' @param track_colors A vector of colors for the annotation tracks. Default is `NULL`.
#' @param track_prop Proportion of the heatmap's height allocated to annotation tracks. Default is 0.1.
#' @param fname Optional. The file name to save the heatmap plot to. Default is `NULL`.
#' @param fwidth Width of the saved heatmap (in inches). Default is 7.
#' @param fheight Height of the saved heatmap (in inches). Default is 5.
#' @param ... Additional arguments to pass to the `ggheatmap` function.
#'
#' @returns A `ggplot` object representing the heatmap.
#' @export
#'
#' @importFrom ggplot2 ggsave
#' @importFrom ggheatmapper ggheatmap add_tracks
#' @importFrom patchwork plot_layout
#'
plt_heatmap <- function(data4hm, annotations = NA,
                        fontsize = 9, track_colors = NULL,
                        track_prop = 0.1,
                        fname = NULL, fwidth = 7, fheight = 5,
                        ...) {
  # Determine the colorbar direction based on the presence of annotations
  colorbar_dir <- if_else(any(is.na(annotations)), "vertical", "horizontal")

  # Generate the basic heatmap using ggheatmap
  hm <- ggheatmapper::ggheatmap(data4hm$table,
                  colv = data4hm$colv,
                  rowv = data4hm$rowv,
                  clustering_method = "ward.D2",
                  fontsize = fontsize,
                  colorbar_dir = colorbar_dir,
                  ...)

  # If annotations are provided, add them to the heatmap
  if (!any(is.na(annotations))) {
    hm <- ggheatmapper::add_tracks(hm,
                     track_columns = annotations,
                     track_colors = track_colors,
                     leg_ncol = 1,
                     track_prop = track_prop,
                     track_pos = "top",
                     legend_action = "collect",
                     fontsize = fontsize) +
      patchwork::plot_layout(guides = "collect")
  }

  # If a file name is provided, save the heatmap as an image
  if (!is.null(fname)) {
    message("-- Saving heatmap at ", fname)
    ggplot2::ggsave(fname, hm, width = fwidth, height = fheight, create.dir = TRUE)
  }

  # Return the heatmap object
  return(hm)
}
