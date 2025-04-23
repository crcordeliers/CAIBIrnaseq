#' Plot Progeny Heatmap
#'
#' Generates a heatmap for Progeny pathway scores, with optional sample annotations and customization options.
#'
#' @param progeny_scores A data frame of pathway scores or a `SummarizedExperiment` object containing pathway scores in its metadata.
#' @param annotations (Optional) Annotations for the samples to be displayed as tracks on the heatmap. Default is `NA`, meaning no annotations are used.
#' @param annotation_prop Proportion of the heatmap width allocated to annotations, if provided. Default is `0.1`.
#' @param annotation_colors (Optional) A named list of colors for the annotation tracks. Names must match the annotation variables.
#' @param fname (Optional) File name to save the heatmap as an image. If `NULL`, the heatmap is displayed but not saved. Default is `NULL`.
#' @param fwidth Width of the saved heatmap image in inches. Default is `7`.
#' @param fheight Height of the saved heatmap image in inches. Default is `5`.
#' @param ... Additional arguments passed to the internal heatmap plotting function.
#'
#' @details
#' If `progeny_scores` is a data frame, it is assumed to contain pathway scores with row names as pathway names and column names as sample IDs.
#' If it is a `SummarizedExperiment` object, pathway scores are extracted from the metadata under `progeny_scores`.
#'
#' Annotations can be added to the heatmap if provided, and their appearance can be customized using `annotation_colors`.
#'
#' @return A heatmap object showing Progeny pathway scores.
#'
#' @importFrom viridisLite viridis
#'
#' @export
plot_progeny_heatmap <- function(progeny_scores,
                                 annotations = NA,
                                 annotation_prop = 0.1,
                                 annotation_colors = NULL,
                                 fname = NULL,
                                 fwidth = 7,
                                 fheight = 5,
                                 ...) {
  if(is.data.frame(progeny_scores)) {
    if(!is.na(annotations)) {
      warning("Plotting pathway score matrix, no sample annotation will be plotted")
    }
    hm_data <- prep_scoredf_hm(progeny_scores)
    hm <- plt_heatmap(hm_data,
                       colors_title = "Pathway score",
                       hm_colors = viridisLite::viridis(100),
                       fname = fname,
                       fwidth = fwidth,
                       fheight = fheight,
                       ...)
  } else {
    exp_data <- progeny_scores
    progeny_scores <- metadata(exp_data)[["progeny_scores"]]
    if(is.null(progeny_scores)) {
      stop('No pathway scores found in `metadata(exp_data)[["progeny_scores"]])`')
    }
    hm_data <- prep_scores_hm(exp_data, progeny_scores)
    hm <- plt_heatmap(hm_data,
                       annotations = annotations,
                       colors_title = "Pathway score",
                       hm_colors = viridisLite::viridis(100),
                       track_prop = annotation_prop,
                       track_colors = annotation_colors,
                       fname = fname,
                       fwidth = 7,
                       fheight = 5,
                       ...)
  }

  return(hm)
}
