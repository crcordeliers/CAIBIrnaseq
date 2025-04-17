#' Plot Heatmap of Progeny Scores with Optional Annotations
#'
#' This function generates a heatmap of progeny scores from either a data frame or a `SummarizedExperiment`
#' object containing progeny scores. Annotations can optionally be added to the heatmap for further customization.
#'
#' @param progeny_scores A `data.frame` or a `SummarizedExperiment` object. If a `data.frame`, it should contain progeny scores. If a `SummarizedExperiment`, it should have progeny scores stored in the metadata.
#' @param annotations Optional. A data frame or list of sample annotations to overlay on the heatmap. Default is `NA` (no annotations).
#' @param annotation_prop Proportion of the heatmap's height to allocate to the annotation tracks. Default is 0.1.
#' @param annotation_colors A vector of colors to use for the annotation tracks. Default is `NULL`.
#' @param fname Optional. The file name to save the heatmap plot to. Default is `NULL` (no file is saved).
#' @param fwidth The width (in inches) of the saved heatmap. Default is 7 inches.
#' @param fheight The height (in inches) of the saved heatmap. Default is 5 inches.
#' @param ... Additional arguments to pass to the `plt_heatmap` function.
#'
#' @returns A `ggplot` object representing the heatmap.
#' @export
#'
#' @importFrom viridisLite viridis
#'
plot_progeny_heatmap <- function(progeny_scores,
                                 annotations = NA,
                                 annotation_prop = 0.1,
                                 annotation_colors = NULL,
                                 fname = NULL,
                                 fwidth = 7,
                                 fheight = 5,
                                 ...) {
  # Si progeny_scores est un data.frame (scores de progeny)
  if (is.data.frame(progeny_scores)) {
    if (!is.na(annotations)) {
      warning("Plotting progeny score matrix, no sample annotation will be plotted")
    }
    # Préparer les données du heatmap pour un data.frame de scores
    hm_data <- prep_scoredf_hm(progeny_scores)
    # Créer et personnaliser le heatmap
    hm <- plt_heatmap(hm_data,
                      colors_title = "Progeny score",
                      hm_colors = viridisLite::viridis(100),
                      fname = fname,
                      fwidth = fwidth,
                      fheight = fheight,
                      ...)
  } else {
    # Si progeny_scores est un SummarizedExperiment, récupérer les scores dans les métadonnées
    exp_data <- progeny_scores
    progeny_scores <- metadata(exp_data)[["progeny_scores"]]

    # Vérifier la présence des scores dans les métadonnées
    if (is.null(progeny_scores)) {
      stop('No progeny scores found in `metadata(exp_data)[["progeny_scores"]])`')
    }

    # Préparer les données du heatmap pour un SummarizedExperiment
    hm_data <- prep_scores_hm(exp_data, progeny_scores)
    # Créer et personnaliser le heatmap
    hm <- plt_heatmap(hm_data,
                      annotations = annotations,
                      colors_title = "Progeny score",
                      hm_colors = viridisLite::viridis(100),
                      track_prop = annotation_prop,
                      track_colors = annotation_colors,
                      fname = fname,
                      fwidth = fwidth,
                      fheight = fheight,
                      ...)
  }

  # Retourner l'objet heatmap
  return(hm)
}
