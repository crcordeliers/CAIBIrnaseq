# CAIBIrnaseq Package Audit

Date: 2026-07-29

Scope: full read-through of all 47 files in `R/`, cross-referenced against the vignettes (`vignettes/*.qmd`) for intended usage. Read-only audit — no fixes applied except where noted as already fixed in a prior session.

## Already fixed (prior session, before this audit)

- `R/pathwayORA.R` — `mutate(padj = p.adjust(pval, ...))` referenced columns (`pval`) that were never created (the data frame only had `PValue`). Renamed ORA's output columns to `pathway`/`pval`/`padj`/`geneRatio`/`bgRatio`/`genes`, matching `pathwayFGSEA()`'s naming.
- `R/plot_pathway_dotplot.R` — ORA-vs-GSEA detection checked for a `geneHits` column that `pathwayORA()` never produces (dead branch). Fixed to check for `geneRatio`.
- `vignettes/venn.qmd` — called `ggsave()` without loading `library(ggplot2)`; added the missing `library()` call.

## Code bugs (certain — will break real usage)

| # | Location | Bug | Impact |
|---|----------|-----|--------|
| 1 | `R/prep_scores_hm.R:54` | Leftover `browser()` call | Every SummarizedExperiment-based pathway/progeny/microenv heatmap hangs in an interactive debugger — breaks any non-interactive render |
| 2 | `R/gexp_preprocess.R:30` | References undefined `gene_mat` (param is `gexp`) | The "auto-lookup annotation via biomaRt" branch is unreachable — errors immediately |
| 3 | `R/gexp_preprocess.R:30` | biomaRt dataset hardcoded to `mmusculus_gene_ensembl` | Human callers silently get mouse gene mappings unless they hand-edit the installed package source |
| 4 | `R/prep_exp_hm.R:44` | References undefined `exp_data` (param is `expData`) | Every `plot_exp_heatmap()` call errors, unless a variable literally named `exp_data` happens to exist in the caller's environment (true in the vignettes — that's why it's never been caught) |
| 5 | `R/pathwayFGSEA.R:35` | `if(is.null(seed)) set.seed(seed)` — inverted logic | Default `seed=0` path never actually seeds (irreproducible FGSEA runs); `seed=NULL` (meant to mean "don't touch RNG") instead calls `set.seed(NULL)`, re-randomizing it |
| 6 | `R/plot_pathway_dotplot.R:31` | `score_name` documented as defaulting to `"resultsORA"` but has no default in the function signature | `plot_pathway_dotplot(exp_data)` throws "argument missing" instead of working as documented |
| 7 | `R/get_annotation_collection.R:49` | Filters on `gs_subcollection`, but msigdbr encodes Hallmark as `gs_collection="H"`/`gs_subcollection=""` (empty, not NA) | `collections = "Hallmark"` — used verbatim in `vignettes/sup_rnaseq.qmd:49` — silently returns **zero** gene sets, no warning. `unsup_rnaseq.qmd` correctly uses `"H"` instead; the two vignettes contradict each other |
| 8 | `R/plot_microenv_heatmap.R:46` | `dir.create(path_dir(fname))` runs unconditionally even though `fname` defaults to `NULL` | Calling with no `fname` (the documented "don't save" default) crashes on `path_dir(NULL)`. Every sibling heatmap function guards this; this one doesn't |
| 9 | `R/zzz.R` (whole file) | `dplyr::desc/relocate/n`, and in `plot_pathway_dotplot.R`: `scale_fill_distiller`, `geom_segment`, `rowwise`, `str_split_fixed` are used but never `importFrom`'d anywhere | A bare `library(CAIBIrnaseq)` session (no `library(tidyverse)` first) breaks `diffExpAnalysis()`, `pathwayFGSEA()`, `pathwayORA()`, `plt_boxplot()`, and `plot_pathway_dotplot()` outright. Only masked because both vignettes attach tidyverse first |
| 10 | `R/zzz.R:29` | `.onLoad` (auto-installs the 4 GitHub-only deps: ggheatmapper, progeny, MCPcounter, mMCPcounter) is entirely commented out | Fresh installs fail the first time any of those deps is needed, with no guidance — worth confirming this was disabled on purpose (e.g. for CRAN) rather than left off by accident |

## Conceptual / statistical concerns

- **`pathwayFGSEA.R:42`** — ranks genes by raw `log2FoldChange` alone, ignoring significance/variance of the estimate. Standard practice ranks by a signed significance statistic (DESeq2's `stat` column, or `-log10(pval)*sign(lfc)`) so a noisy-but-large fold change on a low-count gene doesn't outrank a modest-but-solid one.
- **`score_pathways.R:82`** — docstring says pathways are ranked by "variation," but the code sorts by `rowSums(abs(scores))` (total magnitude), which is a different thing. A pathway that's uniformly strongly active in every sample ranks above one that genuinely swings between conditions.
- **`score_pathways.R:52`** — gene-ID vs. gene-symbol auto-detection needs >100 overlapping IDs with no override. Breaks for anyone using a targeted panel (e.g. a 50-gene NanoString set) even when the annotation is unambiguous.
- **`score_progeny.R`** — unlike `score_pathways()`, doesn't validate gene-ID type before handing the assay to `progeny::progeny()`. If `exp_data` still has Ensembl IDs as rownames (skipping `rebase_gexp` is presented as optional in the vignettes), PROGENy silently matches almost nothing and returns near-empty/NA scores with no warning.

## Heatmap convention

Confirmed: no `pheatmap` usage anywhere in the package or vignettes. `plt_heatmap.R` and all four `prep_*_hm.R` / `plot_*_heatmap.R` wrappers already route through `ggheatmapper` consistently. Nothing to fix here.

## Suggested additional analyses

Given what's already here (DESeq2 DE, ORA/FGSEA, GSVA pathway scoring, PROGENy, MCP-counter/mMCPcounter deconvolution, PCA/clustering), the natural gaps for this pipeline:

- **Batch-effect handling** — nothing currently visualizes or corrects for batch. A thin wrapper around `limma::removeBatchEffect()` (or `sva::ComBat_seq` pre-DESeq2) plus a before/after PCA comparison would slot naturally next to `pca_gexp`/`plot_pca`.
- **Formal sample-outlier detection** — `plot_qc_filters()` exists but there's no explicit outlier flagging (e.g. PCA-distance or sample-correlation based) to catch a bad sample before it contaminates DE/clustering.
- **Alternative DE engine (limma-voom)** — DESeq2 is the only DE method; offering `limma-voom` as an option in `diffExpAnalysis()` would help for very small or unbalanced designs where DESeq2's dispersion estimation is less stable.
- **ORA/FGSEA agreement view** — now that both share `pathway`/`padj` columns, a `plot_pathway_venn()` or similar comparing significant hits between the two methods would be an easy win reusing the existing `plot_venn` machinery.
- **Survival / clinical correlation** — if downstream users have clinical metadata, a simple Kaplan-Meier or Cox wrapper correlating a pathway/gene score against survival would fit naturally alongside `plot_pathway_boxplot`/`plot_pathway_scatter`.
- **Correlation / co-expression heatmap** — a heatmap of top-variable genes or pathways correlated against each other (via `ggheatmapper`, staying on-convention) for quick module/co-regulation exploration.
- **Formal cluster-count selection** — `cluster_k_hc()` exists but nothing tells the user how many clusters to pick. Silhouette score or gap statistic (`cluster::silhouette` / `clusGap`) would pair well with the existing clustering plots.
- **Caching `get_annotation_collection()`** — it re-hits `msigdbr::msigdbr()` (a multi-second download/parse) on every call; memoizing per species would make the pathway-analysis vignette noticeably faster.
- **One-command report** — now that vignettes are `.qmd`, a `render_report(exp_data, ...)` helper that stitches the sup/unsup workflow into a single parameterized Quarto report per dataset would let non-R-fluent collaborators get a full HTML report from one function call.

## Suggested next step

Items 1–9 in "Code bugs" are all small, targeted, high-confidence patches. Item 10 (`.onLoad`) and the conceptual/statistical items involve judgment calls (possibly intentional CRAN-related decisions, or methodological choices) best left to the maintainer.
