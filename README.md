# Multi-Omic Data Integration for Driver Gene Discovery — TCGA-BRCA

## Overview

This project implements a **multi-omic data integration pipeline** to identify **high-confidence driver genes in breast cancer (TCGA-BRCA)** by combining:

* RNA-seq gene expression
* DNA somatic mutation data
* DNA methylation profiles

The workflow integrates R-based preprocessing with Python-based multivariate analysis, using **Canonical Correlation Analysis (CCA)** to detect genes showing coordinated alterations across omic layers.

This pipeline is designed to be **reproducible, modular, and extensible** for other cancer types and multi-omic studies.

---

## Key Findings

The integrated analysis identified two high-confidence driver genes:

### CDH1

* log2 Fold Change: **-2.31**
* Adjusted p-value: **2.5e-51**
* Interpretation: Strong downregulation consistent with **tumor suppressor loss**
* Biological relevance: Known hallmark of **lobular breast carcinoma**

### GATA3

* log2 Fold Change: **+1.55**
* Adjusted p-value: **2.6e-11**
* Interpretation: Significant upregulation
* Biological relevance: Key **luminal subtype transcription factor**

These genes showed **concordant signal across expression, mutation, and methylation**, supporting their role as potential drivers.

---

## Pipeline Architecture

The pipeline combines R (bioconductor ecosystem) and Python (machine learning ecosystem):

| Stage             | Tool                         | Description                                  |
| ----------------- | ---------------------------- | -------------------------------------------- |
| Data Download     | TCGAbiolinks (R)             | Query and download TCGA-BRCA multi-omic data |
| Preprocessing     | SummarizedExperiment (R)     | Filtering, harmonization, imputation         |
| Normalization     | DESeq2, sva (R)              | Variance stabilization and batch correction  |
| Feature Selection | pandas, numpy (Python)       | MAD, variance, mutation recurrence           |
| Integration       | scikit-learn (Python)        | Canonical Correlation Analysis (CCA)         |
| Visualization     | matplotlib, seaborn (Python) | Volcano plot, heatmap, CCA biplot            |

---

## Dataset Summary

| Data Type                          | Features       | Samples |
| ---------------------------------- | -------------- | ------- |
| RNA-seq expression                 | 19,522 genes   | 1,231   |
| DNA mutation                       | 15,359 genes   | 990     |
| DNA methylation                    | 402,482 probes | 895     |
| **Matched samples (intersection)** | —              | **675** |

Only samples with **complete multi-omic coverage** were retained for integration.

---

## Project Structure

```
multiomic-tcga-brca/
│
├── R/
│   ├── 01_download.R
│   ├── 02_preprocess.R
│   └── 03_normalize.R
│
├── python/
│   ├── 04_feature_selection.py
│   ├── 05_cca_integration.py
│   └── 06_visualization.py
│
├── outputs/
│   ├── tables/
│   └── figures/
│
├── data/
│
└── README.md
```

---

## Methodology

### 1. Data Download

* TCGA-BRCA datasets queried using GDC API
* RNA-seq counts
* Somatic mutation MAF files
* DNA methylation (Illumina 450K)

Output stored as SummarizedExperiment objects.

---

### 2. Preprocessing

Steps include:

* Removing low-expression genes
* Mapping Ensembl IDs to gene symbols
* Handling missing values
* Aggregating methylation probes to gene-level
* Mutation binarization (0 = no mutation, 1 = mutated)

---

### 3. Normalization

* RNA-seq normalized using **DESeq2 Variance Stabilizing Transformation (VST)**
* Batch effects corrected using **ComBat (sva)**
* Methylation values standardized
* Mutation matrix left binary

---

### 4. Feature Selection

To reduce dimensionality:

* Median Absolute Deviation (MAD) filtering
* Variance thresholding
* Mutation recurrence filtering
* Top informative genes retained per omic layer

---

### 5. Multi-Omic Integration

Canonical Correlation Analysis (CCA) was applied to:

* Expression matrix
* Mutation matrix
* Methylation matrix

CCA identifies **shared latent components** capturing correlated biological variation across omics.

Driver candidates were selected based on:

* High loading scores
* Consistent directionality
* Statistical significance in differential expression

---

### 6. Visualization

Generated plots include:

* Volcano plot (differential expression)
* Heatmap (top integrated genes)
* CCA biplot (multi-omic relationships)
* Component loadings plot

---

## Installation

### R Dependencies

```
install.packages(c("BiocManager"))
BiocManager::install(c(
  "TCGAbiolinks",
  "DESeq2",
  "sva",
  "SummarizedExperiment",
  "maftools",
  "sesame"
))
```

### Python Dependencies

```
pip install pandas numpy scipy scikit-learn matplotlib seaborn
```

---

## Usage

Run the pipeline sequentially:

### Step 1 — Download Data

```
Rscript R/01_download.R
```

### Step 2 — Preprocess

```
Rscript R/02_preprocess.R
```

### Step 3 — Normalize

```
Rscript R/03_normalize.R
```

### Step 4 — Feature Selection

```
python python/04_feature_selection.py
```

### Step 5 — CCA Integration

```
python python/05_cca_integration.py
```

### Step 6 — Visualization

```
python python/06_visualization.py
```

---

## Output

### Results
![Volcano Plot](outputs/figures/03_volcano.png)
![CCA Biplot](outputs/figures/01_cca_biplot.png)

### Tables

* Driver gene rankings
* CCA component loadings
* Differential expression results

---

## Reproducibility

* Fixed random seeds for CCA
* Modular scripts
* Version-controlled dependencies
* Sample intersection logged

---

## Future Improvements

* Add proteomics data integration
* Replace CCA with MOFA or iCluster
* Pathway enrichment analysis
* Survival analysis using Cox regression
* Docker container for full reproducibility

---

## Citation

If you use this pipeline, please cite:

```
Multi-Omic Integration for Driver Gene Discovery in TCGA-BRCA
Author: Urja Mehta
Year: 2026
```

---

## Author

Urja Mehta
Computer Science | Bioinformatics | Multi-Omics Data Integration

---

## License

MIT License
