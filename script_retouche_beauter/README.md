# 🧬 Subcellular Localization & Isoform Expression Pipeline 🧬

This R pipeline automates the analysis of transcript-level expression and subcellular localization prediction. It integrates:

- 🔹 Transcript TPM values across tissues
- 🔹 Canonical / pseudogene / non-coding classification
- 🔹 Localization predictions (Wolf PSORT)
- 🔹 Experimental validation (Human Protein Atlas)
- 🔹 Functional annotation (nuclear domains via InterProScan)

---

##  Key Features

-  Normalized TPM expression heatmaps
-  Subcellular localization percentage per protein
-  Integration of "approved" / "enhanced" HPA status
-  Automatic detection of nuclear-related domains (zinc finger, helix-turn-helix…)
-  Export of annotated tables as HTML and PDF

---

## Input Files

| 📄 File Name                 | 📌 Description |
|----------------------------|----------------|
| `isoforms_localisation.txt` | Raw Wolf PSORT predictions, format: `details loc:score,...` |
| `hpa_localization.tsv`      | Experimental validation from HPA (`Gene`, `Localization`, `Status`) |
| `results_list_enriched`     | R data frame with expression (TPM), is_canonical, Gene, Transcript, UniProt_ID, etc. |
| `wolf_psort_table_trim27.txt`| Tabulated PSORT output with `Gene|TranscriptID|ProteinID` in column 1 |
| `iprscan5-...tsv`            | InterProScan5 output (TSV, no header) |
| `tsl_lvl.txt`               | Transcript annotation level: `Transcript`, `Raw_TSL`, `GENCODE_raw` |

---

## Output Files

- `sum_isoform.pdf` & `sum_isoform.html`: Expression + annotation summary (TSL, GENCODE basic, etc.)
- `cellular_localization_<gene>.png`: Localization heatmap per gene (prediction + validation)
- `wolf_psort_table_<gene>.pdf`: Table with localization scores and nuclear domain matches

---

## Required R Packages

```r
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(ggplot2)
library(flextable)
library(pagedown)
library(viridis)
