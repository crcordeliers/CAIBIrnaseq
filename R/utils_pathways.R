# Normalize a `pathways` argument into a named list of gene-ID vectors (one
# element per pathway/signature). Accepts either:
#   - a data frame with a `pathway` column and an `id_col` column (e.g. the
#     data frame returned by `get_annotation_collection()`), or
#   - an already-named list of character vectors, e.g.
#     `list(MySignature = c("GENE1", "GENE2"), OtherSignature = c("GENE3"))`,
#     which lets custom gene signatures be scored/tested the same way as an
#     MSigDB collection.
# Shared by pathwayFGSEA(), pathwayORA(), and score_pathways().
.as_gene_sets <- function(pathways, id_col = "gene_symbol") {
  if (is.data.frame(pathways)) {
    if (!all(c("pathway", id_col) %in% colnames(pathways))) {
      stop("`pathways` must be a data frame with a `pathway` column and an `", id_col, "` column ",
           "(or a named list of gene-ID vectors, one per pathway/signature).")
    }
    return(split(pathways[[id_col]], pathways$pathway))
  }
  if (is.list(pathways)) {
    if (is.null(names(pathways)) || any(names(pathways) == "")) {
      stop("`pathways` must be a *named* list of gene-ID vectors (one name per pathway/signature).")
    }
    return(lapply(pathways, as.character))
  }
  stop("`pathways` must be either a data frame or a named list of gene-ID vectors.")
}
