# BioMart and GTEx Data Analysis Pipeline

This repository contains two main pipelines that complement each other for gene expression and transcript analysis:

- **R Pipeline**: Processes gene expression data, calculates expression ratios, generates heatmaps, and prepares summary tables.
- **Python Pipeline**: Processes BioMart data, extracts paralogues, combines transcript lists, filters the GTEx data, and performs hierarchical clustering.

---

## Repository Structure

- **`R_script.R`**: Contains the complete R code for expression analysis, heatmap generation, transcript function analysis, and subcellular localization.
- **`python_script.py`**: Contains the complete Python code for processing BioMart files, filtering GTEx data, and performing hierarchical clustering.
- **`reference/`**: Directory with input files required by both pipelines (e.g., BioMart exports, GTEx data, mapping files).
- **`inter_output/`**: Directory for intermediate outputs from the Python pipeline.
- **Output Files**: Generated PDF/HTML reports, heatmaps, tables, and cluster data files (e.g., `clusterdata_combined.txt`).

---

## Requirements

### R Pipeline Dependencies

- **Packages**: `dplyr`, `tidyr`, `flextable`, `ggplot2`, `pagedown`, `patchwork`, `stringr`, `parallel`, `rlang`, `officer`, `readr`, `png`, `grid`

### Python Pipeline Dependencies

- **Libraries**: `pandas`, `matplotlib`, `seaborn`, `numpy`, `scikit-learn` (for `StandardScaler`), `scipy` (for clustering), plus standard libraries such as `subprocess` and `os`

---
