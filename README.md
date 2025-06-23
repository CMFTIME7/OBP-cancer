# OBP-cancer

A computational and molecular analysis of OBP2A and OBP2B gene and protein expression in cancer, integrating TCGA datasets, mutation profiles, differential gene expression analysis (DGE), and Computational Identification of Biologically Relevant Alterations (CIBRA). This project investigates the functional relevance of odorant-binding proteins in oncology.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data Source: TCGA](https://img.shields.io/badge/Data-TCGA-blue)](https://www.cancer.gov/tcga)
[![R Version](https://img.shields.io/badge/R-4.3.1-blue)](https://cran.r-project.org/)

---

## 🔬 Overview

This repository includes code and data for:
- Expression analysis of OBP2A and OBP2B genes across cancer types (ovarian, breast, uterine, prostate, melanoma, lung, and colorectal cancers).
- Integration of SNV, CNA, and RNA-seq profiles.
- Visualization using `ggplot2`, `maftools`, `EnhancedVolcano`, `ComplexHeatmap`, and more.
- Outputs of oncoplots, boxplots, DGE and CIBRA analysis results.

---

## 📁 Repository Structure

OBP-cancer/
├── data/ # Processed input data from TCGA (SNV, CNA, RNA)
├── scripts/ # Analysis scripts in R
│ ├── visualization/
│ ├── meta-analysis/
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

### 1️⃣ To install R packages:

```r
install.packages(c("tidyverse", "data.table", "ggplot2"))
BiocManager::install(c("clusterProfiler", "fgsea", "org.Hs.eg.db", "ComplexHeatmap"))
``` 

### 2️⃣ Clone the repo:
```r
git clone https://github.com/CMFTIME7/OBP-cancer.git
cd OBP-cancer
```


---
## 🧪 Chemical Sense Linked Figures Overview

This section provides an overview of key figures generated from chemical sense-related data analysis, including CIBRA results, differential gene expression (DGE) analysis, and visual summaries using boxplots and oncoplots.


| 📊 CIBRA                     | 🔬 DGE_analysis                  | 📦 Boxplots                     | 🧬 Oncoplots                    |  
|------------------------------|----------------------------------|---------------------------------|---------------------------------|  
| Overview of binding predictions and receptor interactions using the CIBRA algorithm. | Results of differential gene expression analysis across relevant tissue or cancer datasets. | Visual summaries of expression levels or scores across sample groups. | Mutational landscape highlighting key genes across sample groups. |  
| **Description:** Add insights into the method or figure interpretation. | **Description:** Mention conditions compared, tools used, or notable genes. | **Description:** Add which genes or metrics are being compared. | **Description:** Describe dataset used and how OBPs or related genes appear in cancer profiles. |

---
## Acknowledge
Thanks so much for Soufyan, Mihyeon and Vojta helped with codes and analysis, appreciate Sanne for offer ideas and perspectives as well as feedbacks for project design. And would like to thank Ximeng Wang and Liu Yang for their contributions within their literature overview and M.Sc. research project, respectively. M.C. would like to thank the Giract Flavor research program for a first year PhD bursary (application round 2025). H.M. and M.C. would like to thank the Chinese Scholarship Council for funds (personal PhD grant M.C.) and the Giract Foundation for a 1st year PhD fellowship (personal grant M.C.). H.M. would like to acknowledge the NWO grant. The authors acknowledge the SURFsara compute cluster hosted by SURF and the BAZIS research cluster hosted by VU for the computational time and the provided technical support.
