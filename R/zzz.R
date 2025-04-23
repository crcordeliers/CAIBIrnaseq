utils::globalVariables(c(
  ".",                        # Utilisé dans les pipes (%>%)
  "n",                        # utilisé dans dplyr::summarize
  "gene_id",                  # Utilisé dans rebase_gexp
  "gene_length_kb",           # Utilisé dans rebase_gexp
  "gene_lenght_kb",           # (typo mais conservé de la variable définie par l'utilisateur)
  "gene_description",         # Utilisé dans rebase_gexp
  "gene_biotype",             # Utilisé dans rebase_gexp
  "nfeats",                   # Utilisé dans plot_qc_filters
  "ncounts", "sample_id",     # Utilisé dans plot_qc_filters
  "qc_status",                # Utilisé dans plot_qc_filters
  "gs_collection",            # Utilisé dans get_annotation_collection
  "gs_subcollection",         # Utilisé dans get_annotation_collection
  "gs_name",                  # Utilisé dans get_annotation_collection
  "ensembl_gene",             # Utilisé dans get_annotation_collection
  "gene_symbol",              # Utilisé dans get_annotation_collection
  "collection",               # Utilisé dans get_annotation_collection
  "pathway",                  # Utilisé dans get_annotation_collection
  ".data" ,                   # Utilisé dans get_annotation_collection
  "pathway_scores",           # Utilisé dans cluster_metadata
  "metadata",

  "desc", "log2FoldChange", "PValue", "PAdj", "padj", "Significance",
  "GeneRatio", "GeneRatioNum", "Pathway", "size", "logpadj", "Size",
  "count"
))


# .onLoad <- function(libname, pkgname) {
#   github_pkgs <- list(
#     ggheatmapper = "eurriti/ggheatmapper",
#     progeny = "inmf-lab/progeny",
#     MCPcounter = "ebecht/MCPcounter",
#     mMCPcounter = "diegommcc/mMCPcounter"
#   )
#
#   to_install <- names(github_pkgs)[!vapply(names(github_pkgs), requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
#
#   if (length(to_install) > 0) {
#     if (!requireNamespace("remotes", quietly = TRUE)) {
#       install.packages("remotes")
#     }
#     message("Installing GitHub dependencies for CAIBIrnaseq: ", paste(to_install, collapse = ", "))
#     for (pkg in to_install) {
#       remotes::install_github(github_pkgs[[pkg]])
#     }
#   }
# }
