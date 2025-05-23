<img src="man/figures/CAIBI.png" align="right" alt="" width="120" />

## 🎯 Welcome to CAIBIrnaseq 

**CAIBIrnaseq** is an R package designed to streamline, standardize, and reproduce key steps of gene expression analysis from RNA-seq data. Developed by the **CAIBI platform team**, it integrates powerful tools for:

- **Preprocessing** raw expression data (TPM calculation, filtering, normalization, etc.)
- **Exploratory data analysis** (PCA, clustering, heatmaps)
- **Differential expression analysis** (DESeq2 and visualization of DE genes)
- **Functional and pathway analysis** (FGSEA, ORA, PROGENy scoring)
- **Tumor microenvironment investigation** (MCPcounter, cell-type signatures)


## Quick Installation

In RStudio, try running these command lines :

```
# If not already installed
install.packages("devtools")
install.packages("BiocManager")

# Installation from Github : 
devtools::install_github("crcordeliers/CAIBIrnaseq", dependencies = TRUE)
```

If you have any problems, more informations about installation are available on our package website : 
https://crcordeliers.github.io/CAIBIrnaseq/

**We hope this will help you in your projects !!!**


