utils::globalVariables(c(
  ".",                        # Utilisé dans les pipes (%>%)
    "n",                      # utilisé dans dplyr::summarize
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
  "pathway_scores",            # Utilisé dans cluster_metadata
  "metadata"
))
