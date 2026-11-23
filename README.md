# HEATR3-glioma-eQTL-survival

R code for the study:

> **Germline-regulated glioma risk genes and survival: an integrative GWAS–eQTL–survival study shows that the prognostic signal of HEATR3 in lower-grade glioma is explained by tumour proliferation.**
> Saputri A, Malueka RG, Dwianingsih EK. *[Journal, year, DOI — to be added after publication]*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24468649.svg)](https://doi.org/10.5281/zenodo.24468649)

Zenodo archive: https://doi.org/10.5281/zenodo.24468649

Repository: https://github.com/putriadesagu/HEATR3-glioma-eQTL-survival

---

## What this repository does

The analysis links three layers of public data:

1. **Germline susceptibility.** Genome-wide significant glioblastoma associations from the NHGRI-EBI GWAS Catalog. There were 1,099 records; 64 reached P < 5 × 10⁻⁸, covering 37 lead SNPs and 12 reported genes. Two further candidates were added: *MIR4300HG* (mapped gene of an intergenic variant) and *LMF1* (a non-GBM glioma locus; Melin et al. 2017).
2. **Germline gene regulation.** cis-eQTL evidence in GTEx v8 brain cortex. Six genes were analysed: *HEATR3*, *LMF1*, *EGFR*, *RTEL1*, *SLC16A8* and *MIR4300HG*.
3. **Tumour outcome.** Cox models of overall survival in TCGA-LGG (n = 512), TCGA-GBM (n = 257) and CGGA-693 (grade II–III, n = 240). The models include IDH stratification, time-varying effects, adjustment for a 10-gene proliferation score, and a fixed-effect meta-analysis.

**Main result.** *HEATR3* expression predicts survival in LGG but not GBM, independently of IDH, and the association replicates directionally in CGGA. However, it is fully explained by tumour proliferation in both cohorts: the pooled HR per SD is 1.19 before adjustment for proliferation and 0.99 after.

---

## Repository structure

```
HEATR3-glioma-eQTL-survival/
├── README.md
├── LICENSE
├── scripts/
│   ├── 01_HEATR3_continue.R               # S1 GWAS filtering, S2 GTEx matching, S4 Cox Table (original scale)
│   │                                      #   + HEATR3 x cohort interaction, S5 IDH models and HEATR3 x IDH,
│   │                                      #   S6 Kaplan–Meier (median cut-off, numbers at risk)
│   ├── 02_CGGA_HEATR3_replication.R       # CGGA-693 replication (Table 3), IDH-wildtype subgroup,
│   │                                      #   CGGA time-varying HR, TCGA + CGGA fixed-effect meta-analysis
│   ├── 03_HEATR3_strengthening_project.R  # Proliferation adjustment (TCGA), genome-wide ranking of HEATR3,
│   │                                      #   added prognostic value
│   ├── 04_HEATR3_final_project.R          # Proliferation adjustment in both cohorts + meta (Table 4),
│   │                                      #   time-varying HR before/after proliferation, added value,
│   │                                      #   Ensembl variant identity and LD, Figs 2, 4, 5
│   ├── 05_fix_GBM_Table1_Fig3.R           # Final Table 1 and Fig. 3 per SD (GBM log2(x+1)), ratio of HRs
│   └── 06_make_SuppTables_S1_S2.R         # Supplementary Tables S1, S1b and S2
├── data/                                  # NOT included — see "Input data" below
└── results/                               # created by the scripts
```

Run the scripts in numerical order. Script 04 reads outputs from script 01 (`S4_Cox_Table1.tsv`, S2 audit), and script 06 reads `S1_GWAS_input.tsv` from script 01.

---

## Requirements

Tested with **R 4.2.2** (Windows 10, x86_64).

| Package | Version used | Purpose |
|---|---|---|
| survival | 3.4-0 | Cox models, Schoenfeld tests, time-varying effects, concordance |
| data.table | 1.17.8 | Data handling |
| ggplot2 | 4.0.2 | Figures |
| cowplot | 1.2.0 | Kaplan–Meier panels with risk tables |
| jsonlite | 2.0.0 | Ensembl REST queries (variant identity, LD) |
| fgsea | 1.24.0 | GSEA (optional) |
| msigdbr | — | Gene sets for GSEA (optional; or use MSigDB `.gmt` files) |

```r
install.packages(c("survival", "data.table", "ggplot2", "cowplot", "jsonlite", "msigdbr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("fgsea")
```

An internet connection is needed only for the Ensembl REST step in `02_final_analyses.R`.

---

## Input data (not redistributed)

All data are public. Download them yourself and place them in `data/`, or edit the paths in the `CONFIG` block at the top of each script.

| Data | Source | Files used |
|---|---|---|
| GWAS associations | NHGRI-EBI GWAS Catalog, https://www.ebi.ac.uk/gwas — search "glioblastoma" → Associations → Download | `gwas-association-downloaded_glioblastoma.tsv` (1,099 records) |
| Brain-cortex cis-eQTL | GTEx Portal v8, https://gtexportal.org → Downloads → `GTEx_Analysis_v8_eQTL.tar` | `Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz`, `Brain_Cortex.v8.egenes.txt.gz` |
| TCGA-LGG / TCGA-GBM expression and survival | UCSC Xena, https://xenabrowser.net (TCGA hub) | LGG `HiSeqV2` (log2(RSEM+1)); GBM gene expression; curated survival (OS, OS.time) |
| IDH status (TCGA) | TCGA pan-glioma molecular annotation (Ceccarelli et al. 2016) | merged into the clinical table |
| CGGA-693 | Chinese Glioma Genome Atlas, http://www.cgga.org.cn → Download | `CGGA.mRNAseq_693.RSEM-genes.20200506.txt`, `CGGA.mRNAseq_693_clinical.20200506.txt` |
| Variant identity and LD | Ensembl REST API (GRCh38), https://rest.ensembl.org | queried at run time (1000 Genomes phase 3, EUR) |

**Expected clinical columns** (rename in `CONFIG` → `COLMAP` if yours differ): sample ID, overall-survival time in days, overall-survival event (1 = death), age, sex and, for LGG, IDH status.

---

## How to run

Run the scripts in order from the repository root, either in RStudio or with `Rscript`. Each script writes its tables, figures and a plain-text `SUMMARY*.txt` with every number reported in the manuscript.

Set the project folder in each script (`ROOT` or `BASE` at the top; currently `D:/RESIDENSI/Thesis-PPDS/R`) before running.

```r
source("scripts/01_HEATR3_continue.R")
source("scripts/02_CGGA_HEATR3_replication.R")   # asks for the CGGA files if paths are blank
source("scripts/03_HEATR3_strengthening_project.R")
source("scripts/04_HEATR3_final_project.R")
source("scripts/05_fix_GBM_Table1_Fig3.R")
source("scripts/06_make_SuppTables_S1_S2.R")
```

Steps S7–S10 in script 01 need optional supplementary inputs in `data/` and are not used for the published results; the corresponding analyses were run in scripts 02–04.

---

## Outputs and where they appear in the manuscript

| Output file | Manuscript item |
|---|---|
| `Table1_per_SD_manuscript.csv` | Table 1 — six genes, TCGA-LGG and TCGA-GBM (HR per SD, FDR) |
| `Table2_*` / `P1_time_varying_HR.csv` | Table 2 — HEATR3 IDH-stratified and time-varying models |
| `P1_HEATR3_proliferation_TCGA_CGGA_meta.csv` | Table 3 (CGGA replication, meta-analysis) and Table 4 (proliferation adjustment) |
| `P1_added_value_over_proliferation.csv` | Section 3.5 — C-index and likelihood-ratio test |
| `Fig2_eQTL_pairs.(pdf/tiff)` | Fig. 2 |
| `Fig3_forest_LGG_GBM_perSD.(pdf/tiff)` | Fig. 3 |
| `Fig4_KM_*.(pdf/tiff)` | Fig. 4 (TCGA) and supplementary KM curves (CGGA) |
| `Fig5_forest_proliferation_TCGA_CGGA_meta.(pdf/tiff)` | Fig. 5A |
| `Fig5C_time_varying_HR_TCGA_CGGA.(pdf/tiff)` | Fig. 5B |
| `SuppTable_S1_GWAS_genes.csv`, `SuppTable_S1b_GWAS_associations_all64.csv` | Supplementary Table S1 and file S1b. The script lists the 12 reported genes. In the published table, two rows were added by hand: *MIR4300HG* (mapped gene of intergenic rs11233250) and *LMF1* (non-GBM locus rs3751667, Melin et al. 2017). Implausible odds ratios stored in the Catalog are shown as NR. |
| `P2b_GTEx_variant_identity_verified.csv` / `SuppTable_S2_*.csv` | Supplementary Table S2 |

Figures are saved as vector PDF and 600-dpi LZW-compressed TIFF.

---

## Key analysis choices

- **Expression scale.** TCGA-LGG: log2(RSEM + 1) as in the Xena / merged file. TCGA-GBM: log2(x + 1) for Table 1 and Fig. 3. CGGA-693: RSEM values as provided by CGGA. All are standardised within each cohort, and hazard ratios are per 1 SD.
- **Covariates.** Age and sex. For HEATR3 in LGG, models are stratified by IDH status, because IDH violated proportional hazards.
- **Time-varying effect.** β(t) = β₁ + β₂·[log(1 + t) − log 2], with t in years, so β₁ is the log HR at 1 year. It is tested with a likelihood-ratio test, and HRs are reported at 1, 3 and 5 years.
- **Genome-wide ranking.** The same IDH-stratified, age- and sex-adjusted Cox model is fitted for every expressed TCGA-LGG gene (mean log2 > 1, SD > 0.1), and genes are ranked by |z|.
- **Proliferation score.** Mean z-score of *MKI67, TOP2A, CCNB1, CDK1, BUB1, PCNA, MCM2, AURKA, CCNA2, TPX2*, re-standardised within each cohort.
- **Meta-analysis.** Inverse-variance fixed effect on log HR, with Cochran's Q and I².
- **Multiple testing.** Benjamini–Hochberg FDR within each cohort for the six-gene screen.
- **Kaplan–Meier.** Median split of the whole cohort; log-rank test.

---

## Citation

If you use this code, please cite the article above and the Zenodo archive:

> Saputri A, Malueka RG, Dwianingsih EK. *HEATR3-glioma-eQTL-survival* (v1.0). Zenodo. https://doi.org/10.5281/zenodo.24468649

Please also cite the data resources: GWAS Catalog (Sollis et al. 2023), GTEx Consortium (2020), TCGA, UCSC Xena (Goldman et al. 2020), and CGGA (Zhao et al. 2021).

## License

Code: MIT License (see `LICENSE`). Third-party data remain under the terms of their original providers.

## Contact

Ade Saputri — Department of Anatomical Pathology, Faculty of Medicine, Public Health and Nursing, Universitas Gadjah Mada, Yogyakarta, Indonesia — [ade.saputri@mail.ugm.ac.id]
