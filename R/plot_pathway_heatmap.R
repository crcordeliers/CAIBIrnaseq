#' Plot Pathway Heatmap
#'
#' This function generates a heatmap of pathway scores, with optional sample annotations.
#'
#' @param pathway_scores Either a data frame of pathway scores or a SummarizedExperiment object
#'   containing pathway scores in its metadata.
#' @param pathways A vector specifying the pathways to include in the heatmap. Defaults to the first 20 pathways.
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
#' @importFrom viridisLite viridis
#'
#' @export
#'
plot_pathway_heatmap <- function(pathway_scores,
                                 pathways = 1:20,
                                 annotations = NA,
                                 annotation_prop = 0.1,
                                 annotation_colors = NULL,
                                 fname = NULL,
                                 fwidth = 7,
                                 fheight = 5,
                                 ...) {
  # Si pathway_scores est un data.frame (scores de voies)
  if (is.data.frame(pathway_scores)) {
    if (!is.na(annotations)) {
      warning("Plotting pathway score matrix, no sample annotation will be plotted")
    }
    # Préparer les données du heatmap pour un data.frame de scores
    hm_data <- prep_scoredf_hm(pathway_scores, pathways = pathways)
    # Créer et personnaliser le heatmap
    hm <- plt_heatmap(hm_data,
                      colors_title = "Pathway score",
                      hm_colors = viridisLite::viridis(100),
                      fname = fname,
                      fwidth = fwidth,
                      fheight = fheight,
                      ...)
  } else {
    # Si pathway_scores est un SummarizedExperiment, récupérer les scores dans les métadonnées
    exp_data <- pathway_scores
    pathway_scores <- S4Vectors::metadata(pathway_scores)[["pathway_scores"]]

    # Vérifier la présence des scores dans les métadonnées
    if (is.null(pathway_scores)) {
      stop('No pathway scores found in `metadata(exp_data)[["pathway_scores"]])`')
    }

    # Préparer les données du heatmap pour un SummarizedExperiment
    hm_data <- prep_scores_hm(exp_data, pathway_scores, pathways)
    # Créer et personnaliser le heatmap
    hm <- plt_heatmap(hm_data,
                      annotations = annotations,
                      colors_title = "Pathway score",
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
