# Latent Genomic Factors and and Treatment Outcomes

This repository contains scripts for preprocessing GWAS summary statistics, computing polygenic risk scores (PRS) using PRS-CS, preparing cohort-level datasets with PRS and covariates, and testing their association with treatment outcomes.

The code is designed to be applied across different cohorts and can be adapted to datasets with similar structure.

The pipeline integrates:
- GWAS summary statistics preprocessing  
- PRS computation (PRS-CS + PLINK2)  
- Dataset preparation  
- Logistic regression analyses  
- Multiple testing correction  

---

## Repository structure

```
scripts_PRSs_latent_factors_treatment_outcomes/
│
├── 00_compute_effective_sample_size.sh
├── 01_munge_sumstats.sh
├── 02_run_PRScs.sh
├── 03_compute_PRSs.sh
├── 04_prepare_dataset.R
├── 05_run_logistic_models.R
├── 06_FDR_correction_multiple_testing.R
```

---

## Pipeline overview

The analysis is organised into sequential steps:

### 1. Summary statistics preprocessing
- `00_compute_effective_sample_size.sh`  
  Computes the effective sample size (Neff) from GWAS summary statistics.

- `01_munge_sumstats.sh`  
  Formats and filters summary statistics using LDSC (e.g., INFO and MAF thresholds, HapMap3 SNP restriction).

---

### 2. PRS computation
- `02_run_PRScs.sh`  
  Runs PRS-CS separately across chromosomes using a Bayesian framework and parallel tmux sessions.

- `03_compute_PRSs.sh`  
  Combines chromosome-specific weights and computes PRS using PLINK2.

---

### 3. Dataset preparation
- `04_prepare_dataset.R`  
  - Reads a cohort-level dataset  
  - Merges PRS with the cohort dataset  
  - Merges external covariates (e.g., principal components, site variables)  
  - Harmonises duplicated variables created during merging  
  - Converts variable types where needed  
  - Computes z-scores for PRS  
  - Derives analysis-ready outcomes when required  

---

### 4. Association analyses
- `05_run_logistic_models.R`  
  Runs logistic regression models across multiple outcomes and PRS, including interaction terms with selected moderators (e.g., age, sex, baseline severity).

Outputs include:
- Model summaries  
- Odds ratios and confidence intervals  
- Log files  
- Plots for PRS-related effects meeting the specified p-value threshold  

---

### 5. Multiple testing correction
- `06_FDR_correction_multiple_testing.R`  
  Adjusts p-values using the Benjamini–Hochberg false discovery rate.

Multiple testing is controlled across outcome-specific tests for each PRS.  
Significance is defined as:

```
q < 0.10 (two-sided)
```

---

## Requirements

### Software
- Bash (Unix-based system recommended)
- R (≥ 4.0)
- Python 2 (for LDSC)
- Python 3 (for PRS-CS)
- PLINK2
- tmux

### External tools
- LDSC: https://github.com/bulik/ldsc  
- PRS-CS: https://github.com/getian107/PRScs  

### R packages
- dplyr  
- readr  
- stringr  
- ggplot2  
- broom  

---

## Data

- GWAS summary statistics are publicly available and should be downloaded separately.
- Individual-level cohort data are not included in this repository.

All scripts require user-defined paths to input and output files.

---

## Usage

Run the pipeline in the following order:

```
00_compute_effective_sample_size.sh
01_munge_sumstats.sh
02_run_PRScs.sh
03_compute_PRSs.sh
04_prepare_dataset.R
05_run_logistic_models.R
06_FDR_correction_multiple_testing.R
```

---

## Notes

- File paths must be specified manually in each script.  
- Scripts assume consistent subject identifiers across datasets (e.g., genetic ID and PRS ID) after harmonisation if needed.  
- PRS are standardised (z-scores) before analysis.  
- Logistic models are fitted using complete-case data.  
- Variable names, outcome definitions, and covariate structures may need to be adapted depending on the cohort.  

---

## Contact

For questions or collaboration, please contact the repository owner.
