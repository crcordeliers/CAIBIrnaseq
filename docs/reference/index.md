# Package index

## Preprocessing

- [`filter_gexp()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/filter_gexp.md)
  : Filter Gene Expression Data
- [`normalize_gexp()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/normalize_gexp.md)
  : Normalize Gene Expression Data
- [`gexp_preprocess()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/gexp_preprocess.md)
  : Pre-process a gene expression matrix to have different gene
  annotation
- [`rebase_gexp()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/rebase_gexp.md)
  : Rebase Gene Expression Matrix
- [`fixMarsExcel()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/fixMarsExcel.md)
  : Fix Excel-Corrupted Gene Names
- [`robust_cv()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/robust_cv.md)
  : Robust Coefficient of Variation
- [`highly_variable_genes()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/highly_variable_genes.md)
  : Identify Highly Variable Genes
- [`estimateTPM()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/estimateTPM.md)
  : Estimate Transcripts Per Million (TPM)
- [`readsPerMillion()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/readsPerMillion.md)
  : Normalize Counts to Reads Per Million (RPM)
- [`transcripts_per_million()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/transcripts_per_million.md)
  : Transcripts Per Million (TPM) Calculation

## Exploratory Analysis

- [`pca_gexp()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/pca_gexp.md)
  : Perform PCA on Gene Expression Data
- [`plot_pca()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_pca.md)
  : Plot PCA Results for Gene Expression Data
- [`cluster_exp()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/cluster_exp.md)
  : Cluster Samples Based on Gene Expression
- [`cluster_k_hc()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/cluster_k_hc.md)
  : Hierarchical Clustering with Optional PCA Dimensionality Reduction
- [`cluster_metadata()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/cluster_metadata.md)
  : Cluster Metadata in Expression Data

## Differential Expression Analysis

- [`diffExpAnalysis()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/diffExpAnalysis.md)
  : Differential expression analysis (DESeq2 / limma-voom wrapper)
- [`plot_exp_volcano()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_exp_volcano.md)
  : Plot Volcano Plot of Differential Expression Results
- [`plot_exp_boxplot()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_exp_boxplot.md)
  : Plot Expression Boxplot (with validation)
- [`plot_exp_heatmap()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_exp_heatmap.md)
  : Plot Expression Heatmap
- [`plot_exp_scatter()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_exp_scatter.md)
  : Plot Expression Scatter Plot
- [`plot_qc_filters()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_qc_filters.md)
  : Plot QC Filters for Gene Expression Data
- [`plot_venn()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_venn.md)
  : Plot a Venn diagram for two or three vectors

## Pathway Analysis

- [`pathwayAnalysis()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/pathwayAnalysis.md)
  : Pathway Analysis
- [`pathwayFGSEA()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/pathwayFGSEA.md)
  : Pathway Enrichment Using FGSEA
- [`pathwayORA()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/pathwayORA.md)
  : Pathway Over-Representation Analysis (ORA)
- [`score_pathways()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/score_pathways.md)
  : Score Pathways
- [`score_progeny()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/score_progeny.md)
  : Score PROGENy Pathways
- [`get_pathway_df()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/get_pathway_df.md)
  : Extract pathway scores and metadata from a SummarizedExperiment
- [`get_pathway_genes()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/get_pathway_genes.md)
  : Get Pathway Genes
- [`plot_pathway_dotplot()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_pathway_dotplot.md)
  : Plot Pathway Analysis Results as a Dot Plot
- [`plot_pathway_heatmap()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_pathway_heatmap.md)
  : Plot Pathway Heatmap
- [`plot_pathway_boxplot()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_pathway_boxplot.md)
  : Plot Pathway Boxplot
- [`plot_pathway_scatter()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_pathway_scatter.md)
  : Plot Pathway Scatter Plot

## Microenvironment & Cell-type Scores

- [`mcp_counter()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/mcp_counter.md)
  : MCP-counter wrapper for SummarizedExperiment objects
- [`plot_microenv_heatmap()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_microenv_heatmap.md)
  : Plot Microenvironment Heatmap
- [`plot_progeny_heatmap()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plot_progeny_heatmap.md)
  : Plot Progeny Heatmap

## Utilities and Helpers

- [`get_exp_df()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/get_exp_df.md)
  : Extract gene expression and sample annotations
- [`get_annotation_collection()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/get_annotation_collection.md)
  : Get Annotation Collection
- [`getEmptyRows()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/getEmptyRows.md)
  : Get Empty Rows in a Matrix
- [`plt_boxplot()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plt_boxplot.md)
  : Boxplot of gene expression by annotation group
- [`plt_heatmap()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plt_heatmap.md)
  : Plot Heatmap with Optional Annotations and Customization
- [`plt_scatter()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/plt_scatter.md)
  : Scatter plot of gene expression
- [`prep_exp_hm()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/prep_exp_hm.md)
  : Prepare Gene Expression Data for Heatmap Visualization
- [`prep_scoredf_hm()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/prep_scoredf_hm.md)
  : Prepare Pathway Scores Data Frame for Heatmap
- [`prep_scores_hm()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/prep_scores_hm.md)
  : Prepare Pathway Scores for Heatmap Visualization
- [`read_rnaseq_out()`](https://crcordeliers.github.io/CAIBIrnaseq/reference/read_rnaseq_out.md)
  : Read RNA-seq Output from nf-core/rnaseq
