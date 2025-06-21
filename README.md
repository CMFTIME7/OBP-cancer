# OBP-cancer

A computational and molecular analysis of OBP2A and OBP2B gene expression in cancer, integrating TCGA datasets, mutation profiles, and gene set enrichment analysis (GSEA). This project investigates the functional relevance of odorant-binding proteins in oncology.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data Source: TCGA](https://img.shields.io/badge/Data-TCGA-blue)](https://www.cancer.gov/tcga)
[![R Version](https://img.shields.io/badge/R-4.3.1-blue)](https://cran.r-project.org/)

## Acknowledge
Thanks so much for Soufyan, Mihyeon and Vojta helped with codes and analysis, appreciate Sanne for offer ideas and perspectives as well as feedbacks for project design. And would like to thank Ximeng Wang and Liu Yang for their contributions within their literature overview and M.Sc. research project, respectively. M.C. would like to thank the Giract Flavor research program for a first year PhD bursary (application round 2025). H.M. and M.C. would like to thank the Chinese Scholarship Council for funds (personal PhD grant M.C.) and the Giract Foundation for a 1st year PhD fellowship (personal grant M.C.). H.M. would like to acknowledge the NWO grant. The authors acknowledge the SURFsara compute cluster hosted by SURF and the BAZIS research cluster hosted by VU for the computational time and the provided technical support.

---

## 🔬 Overview

This repository includes code and data for:
- Expression analysis of OBP2A and OBP2B genes across cancer types.
- Integration of SNV, CNA, and RNA-seq profiles.
- Enrichment analysis using `clusterProfiler`, `fgsea`, and MSigDB.
- Visualization using `ggplot2`, `ComplexHeatmap`, and more.

---

## 📁 Repository Structure

OBP-cancer/
├── data/ # Processed input data from TCGA (SNV, CNA, RNA)
├── scripts/ # Analysis scripts in R
│ ├── preprocessing/
│ ├── gsea/
│ ├── visualization/
├── results/ # Output plots, tables, and figures
├── notebooks/ # Optional Jupyter/Rmarkdown exploratory analysis
├── docs/ # Documentation (setup, methods)
├── .github/ # GitHub-specific templates
│ └── ISSUE_TEMPLATE.md
├── LICENSE
├── README.md
└── CITATION.cff


---

## 🧰 Requirements

- R ≥ 4.3.1
- R packages: `tidyverse`, `clusterProfiler`, `fgsea`, `org.Hs.eg.db`, `ggplot2`, `data.table`, `ComplexHeatmap`
- Optional: `Python ≥ 3.10` for cross-platform processing or plotting

To install R packages:

```r
install.packages(c("tidyverse", "data.table", "ggplot2"))
BiocManager::install(c("clusterProfiler", "fgsea", "org.Hs.eg.db", "ComplexHeatmap"))


# Clone the repo
git clone https://github.com/yourusername/OBP-cancer.git
cd OBP-cancer

# Run preprocessing
Rscript scripts/preprocessing/load_and_merge_data.R

# Run GSEA
Rscript scripts/gsea/run_msigdb_enrichment.R

# Generate plots
Rscript scripts/visualization/plot_heatmaps.R
---
## 🧪 Chemical Sense Linked Figures Overview

This section provides an overview of key figures generated from chemical sense-related data analysis, including CIBRA results, differential gene expression (DGE) analysis, and visual summaries using boxplots and oncoplots.

### 📊 CIBRA  
Overview of binding predictions and receptor interactions using the CIBRA algorithm.  
_Description: Add insights into the method or figure interpretation._

### 🔬 DGE_analysis  
Results of differential gene expression analysis across relevant tissue or cancer datasets.  
_Description: Mention conditions compared, tools used, or notable genes._

### 📦 Boxplots  
Visual summaries of expression levels or scores across sample groups.  
_Description: Add which genes or metrics are being compared._

### 🧬 Oncoplots  
Mutational landscape highlighting key genes across sample groups.  
_Description: Describe dataset used and how OBPs or related genes appear in cancer profiles._
