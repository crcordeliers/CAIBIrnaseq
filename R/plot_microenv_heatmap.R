#' Plot Microenvironment Heatmap
#'
#' This function generates a heatmap of microenvironment population scores, with optional sample annotations.
#'
#' @param exp_data A SummarizedExperiment object containing microenvironment population scores in its metadata.
#' @param annotations Optional. Sample annotations to include in the heatmap. Defaults to NA.
#' @param annotation_prop A numeric value specifying the proportion of the heatmap allocated to annotations. Defaults to 0.1.
#' @param annotation_colors Optional. A list of colors for annotations. Defaults to NULL.
#' @param fname Optional. A character string specifying the file name to save the heatmap. Defaults to NULL.
#' @param fwidth Numeric. The width of the saved heatmap file. Defaults to 7.
#' @param fheight Numeric. The height of the saved heatmap file. Defaults to 5.
#' @param ... Additional parameters passed to the heatmap plotting function.
#'
#' @return A heatmap plot object.
#'
#' @importFrom viridisLite plasma
#'
#' @export
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
  microenv_scores <- S4Vectors::metadata(exp_data)[["microenv_scores"]]

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
