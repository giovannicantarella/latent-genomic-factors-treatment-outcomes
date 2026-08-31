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

```
scripts_PRSs_latent_factors_treatment_outcomes/
├── 00_compute_effective_sample_size.sh
├── 01_munge_sumstats.sh
├── 02_run_PRScs.sh
├── 03_compute_PRSs.sh
├── 04_prepare_dataset.R
├── 05_run_frequentist_models.R
├── 06_FDR_correction_multiple_testing.R
├── 07_run_bayesian_sensitivity_analyses.R
```

Each script represents one step of the pipeline and can be adapted to different cohorts by modifying input data and variable mappings.

---

## Analytical workflow

The analysis is structured in eight sequential steps:

1. Effective sample size calculation  
2. GWAS summary statistics preprocessing  
3. PRS-CS estimation  
4. Individual-level PGS computation  
5. Dataset preparation  
6. Frequentist association analyses  
7. Multiple testing correction  
8. Bayesian sensitivity analyses  

---

## Scripts description

### 00 – Effective sample size calculation

Computes the effective sample size (Neff) from GWAS summary statistics.

Main steps:
- Read GWAS summary statistics  
- Calculate effective sample size  
- Generate output for downstream preprocessing  

---

### 01 – Summary statistics preprocessing

Formats and filters summary statistics using LDSC.

Main features:
- INFO filtering  
- MAF filtering  
- HapMap3 SNP restriction  
- LDSC-compatible formatting  

Output:
- Filtered summary statistics for downstream PGS computation  

---

### 02 – PRS-CS

Runs PRS-CS separately across chromosomes.

Main features:
- Bayesian continuous-shrinkage framework  
- Chromosome-wise estimation  
- Parallel execution using tmux sessions  
- Posterior SNP weight estimation  

Output:
- Chromosome-specific posterior SNP weights  

---

### 03 – PGS computation

Combines chromosome-specific weights and computes individual-level PGSs using PLINK2.

Main steps:
- Combine chromosome-specific PRS-CS weights  
- Prepare scoring files  
- Compute individual-level PGSs  
- Generate score files for each latent genomic factor  

---

### 04 – Dataset preparation

Builds the cohort-level analytic dataset.

Main steps:
- Read cohort-level data  
- Merge PGSs with the cohort dataset  
- Merge external covariates (e.g. ancestry principal components and site variables)  
- Harmonise duplicated variables created during merging  
- Convert variable types where needed  
- Standardise PGSs  
- Derive analysis-ready outcomes when required  

---

### 05 – Frequentist association analyses

Runs logistic regression models across multiple outcomes and PGSs.

Main features:
- Base models including PGS and prespecified covariates  
- PGS × baseline severity interaction models  
- PGS × sex interaction models  
- PGS × age interaction models  
- Model comparison using likelihood-ratio tests and AIC  
- Interaction retention based on LRT p < 0.05 and ΔAIC ≥ 2  
- 10-fold out-of-fold AUC estimation  
- Post-hoc predicted probabilities for retained interactions  
- Simple PGS slopes and contrasts across moderator levels  

Output:
- Model selection summaries  
- Odds ratios and 95% confidence intervals  
- Out-of-fold AUC estimates  
- Predicted probabilities  
- Simple-slope estimates and contrasts  

---

### 06 – Multiple testing correction

Applies the Benjamini–Hochberg false discovery rate correction to frequentist association results.

Main features:
- Correction across outcome-specific tests within each PGS  
- Two-sided statistical testing  
- FDR significance threshold of q < 0.05  

---

### 07 – Bayesian sensitivity analyses

Re-estimates the models selected in the frequentist analysis using Bayesian logistic regression.

Main features:
- Same model specification selected in the frequentist analysis  
- No additional model-selection step  
- Weakly informative Normal(0, 0.5) priors for regression coefficients  
- Posterior odds ratios and 95% credible intervals  
- Posterior predicted probabilities for retained interaction models  

Output:
- Posterior coefficient estimates  
- Posterior odds ratios and credible intervals  
- Posterior predicted probabilities  

---

## Requirements

R (≥ 4.0 recommended)

Additional software:
- Bash  
- Python 2 (for LDSC)  
- Python 3 (for PRS-CS)  
- PLINK2  
- tmux  

Key R packages:
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

External tools:
- LDSC: https://github.com/bulik/ldsc  
- PRS-CS: https://github.com/getian107/PRScs  

---

## Data

- GWAS summary statistics are publicly available and should be downloaded separately  
- Individual-level cohort data are not included in this repository  
- Input and output paths must be adapted locally  

---

## General principles

- Scripts are generic templates, not cohort-specific pipelines  
- Cohort differences are handled at the level of:
  - input data  
  - variable naming  
  - outcome definitions  
  - covariate and site structures  
- The same complete-case sample is used when comparing candidate models for each PGS-outcome pair  
- Continuous predictors are standardised before regression modelling  
- Bayesian models reproduce the specifications selected in the frequentist analysis  

---

## Notes

- Scripts are intended to be run sequentially  
- File paths must be adapted locally  
- Subject identifiers should be harmonised across datasets before merging  
- Raw individual-level data are not included in this repository  
- Outputs are generated as tables and model objects for downstream use  

---

## Contact

For questions or collaboration, please contact the repository owner.
````
