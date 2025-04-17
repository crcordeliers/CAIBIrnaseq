#' Plot Heatmap of Microenvironment Population Scores with Optional Annotations
#'
#' This function generates a heatmap of microenvironment population scores from a `SummarizedExperiment`
#' object, with optional annotations that can be added to the heatmap for additional customization.
#'
#' @param exp_data A `SummarizedExperiment` object containing microenvironment population scores stored in the metadata.
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
#' @importFrom viridisLite plasma
#'
plot_microenv_heatmap <- function(exp_data,
                                  annotations = NA,
                                  annotation_prop = 0.1,
                                  annotation_colors = NULL,
                                  fname = NULL,
                                  fwidth = 7,
                                  fheight = 5,
                                  ...) {
  # Extraire les scores de microenvironnement depuis les métadonnées
  microenv_scores <- metadata(exp_data)[["microenv_scores"]]

  # Vérifier si les scores existent dans les métadonnées
  if (is.null(microenv_scores)) {
    stop('No microenvironment population scores found in `metadata(exp_data)[["microenv_scores"]])`')
  }

  # Préparer les données pour le heatmap en utilisant la fonction `prep_scores_hm`
  hm_data <- prep_scores_hm(exp_data, microenv_scores)

  # Générer le heatmap en utilisant `plt_heatmap`
  hm <- plt_heatmap(hm_data,
                    center = TRUE,  # Centrer les données pour le heatmap
                    scale = TRUE,   # Normaliser les données
                    annotations = annotations,
                    colors_title = "Population score",  # Titre des couleurs
                    hm_colors = viridisLite::plasma(100),  # Choisir une palette de couleurs
                    track_prop = annotation_prop,  # Proportion des annotations
                    track_colors = annotation_colors,  # Couleurs des annotations
                    fname = fname,  # Nom du fichier pour enregistrer l'image
                    fwidth = fwidth,  # Largeur de l'image
                    fheight = fheight,  # Hauteur de l'image
                    ...)  # Passer d'autres arguments à `plt_heatmap`

  # Retourner l'objet heatmap
  return(hm)
}
