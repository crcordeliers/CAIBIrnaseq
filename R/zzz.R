# utils::globalVariables(c(
#   ".",                        # Used in pipes (%>%)
#   "n",                        # Used in dplyr::summarize
#   "gene_id",                  # Used in rebase_gexp
#   "gene_length_kb",           # Used in rebase_gexp
#   "gene_lenght_kb",           # (typo but kept from user-defined variable)
#   "gene_description",         # Used in rebase_gexp
#   "gene_biotype",             # Used in rebase_gexp
#   "nfeats",                   # Used in plot_qc_filters
#   "ncounts", "sample_id",     # Used in plot_qc_filters
#   "qc_status",                # Used in plot_qc_filters
#   "gs_collection",            # Used in get_annotation_collection
#   "gs_subcollection",         # Used in get_annotation_collection
#   "gs_name",                  # Used in get_annotation_collection
#   "ensembl_gene",             # Used in get_annotation_collection
#   "gene_symbol",              # Used in get_annotation_collection
#   "collection",               # Used in get_annotation_collection
#   "pathway",                  # Used in get_annotation_collection
#   ".data" ,                   # Used in get_annotation_collection
#   "pathway_scores",           # Used in cluster_metadata
#   "metadata",
#
#   "desc", "log2FoldChange", "PValue", "PAdj", "padj", "Significance",
#   "GeneRatio", "GeneRatioNum", "Pathway", "size", "logpadj", "Size",
#   "count", ":="
# ))


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
