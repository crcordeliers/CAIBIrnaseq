library(tidyverse)
library(ggrepel)
deseq_volcanoPlot <- function(resLFC,
                              contrast_name = "",
                              lfc_tresh = 2.5,
                              padj_tresh = 0.05,
                              ngenes_show = 10,
                              gene_col = "external_gene_name",
                              annot_pos = c(Inf, Inf)) {
  #-- Get annotations
  #---- Genes to show
  genes_up <- resLFC %>%
    filter(log2FoldChange > lfc_tresh) %>%
    slice_min(order_by = padj, n = ngenes_show/2) %>%
    pull({{gene_col}})
  genes_down <- resLFC %>%
    filter(log2FoldChange < -lfc_tresh) %>%
    slice_min(order_by = padj, n = ngenes_show/2) %>%
    pull({{gene_col}})
  genes_label <- c(genes_up, genes_down)
  #---- x-axis
  fc_leg <- paste("FC", contrast_name)
  #---- Number of DEGs
  ndegs <- sum(resLFC$padj < padj_tresh, na.rm = TRUE)
  deg_annot <-  paste0("atop(bold('DEGs:'), ", ndegs, "/", nrow(resLFC), ")")

  #-- Plot data
  data2plot <- resLFC %>%
    mutate(lpval = -log10(padj),
           signif = case_when(
             log2FoldChange > lfc_tresh & padj < padj_tresh ~ "up",
             log2FoldChange < -lfc_tresh & padj < padj_tresh ~ "down",
             padj < 0.05 ~ "signif",
             TRUE ~ "ns"
           ),
           gene_name = pull(resLFC, {{gene_col}}),
           label = ifelse(gene_name %in% genes_label, gene_name, NA))

  #-- Get annotation position params
  hjust_value <- case_when(
    annot_pos[1] == Inf ~ 1,
    annot_pos[1] == -Inf ~ -1,
    TRUE ~ 0
  )
  vjust_value <- case_when(
    annot_pos[2] == Inf ~ 1,
    annot_pos[2] == -Inf ~ -1,
    TRUE ~ 0
  )

  plt <- ggplot(data2plot, aes(log2FoldChange, lpval)) +
    geom_point(aes(color = signif), alpha = 0.7) +
    geom_text_repel(aes(label = label), size = 3) +
    annotate("text", x=annot_pos[1], y=annot_pos[2],
             vjust=vjust_value, hjust=hjust_value, label=deg_annot,
             size = 3.5, parse = TRUE) +
    scale_color_manual(values = c("up" = "#ef8a62",
                                  "down" = "#67a9cf",
                                  "signif" = "grey60",
                                  "ns" = "grey90")) +
    labs(x = bquote("log"[2]*.(fc_leg)), y = expression(-"log"[10]*"adj. p-value")) +
    guides(color = "none") +
    theme_minimal()

  return(plt)

}
