````markdown
# Latent Genomic Factors and Treatment Outcomes

This repository contains scripts for preprocessing GWAS summary statistics, computing polygenic scores (PGSs) using PRS-CS, preparing cohort-level datasets with PGSs and covariates, and testing their association with treatment outcomes.

The code is designed to be applied across different cohorts and can be adapted to datasets with similar structure.

The pipeline integrates:
- GWAS summary statistics preprocessing  
- PGS computation (PRS-CS + PLINK2)  
- Dataset preparation  
- Frequentist logistic regression analyses  
- Multiple testing correction  
- Bayesian sensitivity analyses  

---

## Repository structure

```text
scripts_PRSs_latent_factors_treatment_outcomes/
│
├── 00_compute_effective_sample_size.sh
├── 01_munge_sumstats.sh
├── 02_run_PRScs.sh
├── 03_compute_PRSs.sh
├── 04_prepare_dataset.R
├── 05_run_frequentist_models.R
├── 06_FDR_correction_multiple_testing.R
└── 07_run_bayesian_sensitivity_analyses.R
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

### 2. PGS computation
- `02_run_PRScs.sh`  
  Runs PRS-CS separately across chromosomes using a Bayesian framework and parallel tmux sessions.

- `03_compute_PRSs.sh`  
  Combines chromosome-specific weights and computes PGSs using PLINK2.

---

### 3. Dataset preparation
- `04_prepare_dataset.R`  
  - Reads a cohort-level dataset  
  - Merges PGSs with the cohort dataset  
  - Merges external covariates (e.g., principal components, site variables)  
  - Harmonises duplicated variables created during merging  
  - Converts variable types where needed  
  - Computes standardised PGSs  
  - Derives analysis-ready outcomes when required  

---

### 4. Frequentist association analyses
- `05_run_frequentist_models.R`  
  Runs logistic regression models across multiple outcomes and PGSs, including interaction terms with selected moderators (e.g., age, sex, baseline severity).

Candidate interaction models are compared with the base model using likelihood-ratio tests and AIC.

Interaction terms are retained when both criteria are met:
- LRT `p < 0.05`  
- `ΔAIC ≥ 2`  

Outputs include:
- Model selection summaries  
- Odds ratios and confidence intervals  
- Out-of-fold AUC estimates  
- Predicted probabilities for retained interactions  
- Simple-slope analyses  
- Contrasts between PGS slopes across moderator levels  

---

### 5. Multiple testing correction
- `06_FDR_correction_multiple_testing.R`  
  Adjusts p-values using the Benjamini–Hochberg false discovery rate.

Multiple testing is controlled across outcome-specific tests for each PGS.  
Significance is defined as:

```text
q < 0.05 (two-sided)
```

---

### 6. Bayesian sensitivity analyses
- `07_run_bayesian_sensitivity_analyses.R`  
  Re-estimates the models selected in the frequentist analysis using Bayesian logistic regression.

Weakly informative Normal(0, 0.5) priors are applied to regression coefficients on the log-odds scale.

Outputs include:
- Posterior odds ratios  
- 95% credible intervals  
- Posterior predicted probabilities for retained interaction models  

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
- tidyr  
- tibble  
- readr  
- stringr  
- ggplot2  
- emmeans  
- writexl  
- rstanarm  
- posterior  

---

## Data

- GWAS summary statistics are publicly available and should be downloaded separately.
- Individual-level cohort data are not included in this repository.

All scripts require user-defined paths to input and output files.

---

## Usage

Run the pipeline in the following order:

```text
00_compute_effective_sample_size.sh
01_munge_sumstats.sh
02_run_PRScs.sh
03_compute_PRSs.sh
04_prepare_dataset.R
05_run_frequentist_models.R
06_FDR_correction_multiple_testing.R
07_run_bayesian_sensitivity_analyses.R
```

---

## Notes

- File paths must be specified manually in each script.  
- Scripts assume consistent subject identifiers across datasets after harmonisation if needed.  
- Continuous predictors are standardised before analysis.  
- Logistic models are fitted using complete-case data.  
- The same complete-case sample is used when comparing candidate models for each PGS-outcome pair.  
- Variable names, outcome definitions, and covariate structures may need to be adapted depending on the cohort.  
- Bayesian models reproduce the model specification selected in the frequentist analysis.

---

## Contact

For questions or collaboration, please contact the repository owner.
````
