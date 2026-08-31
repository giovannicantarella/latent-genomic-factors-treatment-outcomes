````markdown
# Latent Genomic Factors and Treatment Outcomes

This repository contains scripts for preprocessing GWAS summary statistics, computing polygenic scores (PGSs) using PRS-CS, preparing cohort-level datasets with PGSs and covariates, and testing their associations with treatment outcomes.

The code is designed to be applied across different cohorts and can be adapted to datasets with a similar structure.

The pipeline includes:
- GWAS summary statistics preprocessing  
- PGS computation (PRS-CS + PLINK2)  
- Dataset preparation  
- Frequentist logistic regression analyses  
- Model comparison and selection  
- Out-of-fold model discrimination  
- Post-hoc interaction analyses  
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

The analysis is organised into sequential steps.

### 1. Summary statistics preprocessing

- `00_compute_effective_sample_size.sh`  
  Computes the effective sample size (Neff) from GWAS summary statistics.

- `01_munge_sumstats.sh`  
  Formats and filters GWAS summary statistics using LDSC, including INFO and MAF filtering and restriction to HapMap3 SNPs.

---

### 2. PGS computation

- `02_run_PRScs.sh`  
  Runs PRS-CS separately across chromosomes using a Bayesian continuous-shrinkage framework and parallel tmux sessions.

- `03_compute_PRSs.sh`  
  Combines chromosome-specific posterior SNP weights and computes individual-level PGSs using PLINK2.

---

### 3. Dataset preparation

- `04_prepare_dataset.R`  
  - Reads a cohort-level clinical dataset  
  - Reads and merges PGS files  
  - Merges external covariates, such as ancestry principal components and site variables  
  - Harmonises duplicated variables created during merging  
  - Converts variable types where needed  
  - Standardises PGSs  
  - Derives analysis-ready outcomes when required  

---

### 4. Frequentist association analyses

- `05_run_frequentist_models.R`  

  Fits multivariable logistic regression models across multiple binary outcomes and PGSs.

For each PGS-outcome pair, the script fits:
- A base model containing the PGS and prespecified covariates  
- A PGS × baseline severity interaction model  
- A PGS × sex interaction model  
- A PGS × age interaction model  

All models include the prespecified clinical, demographic, ancestry, and site covariates.

Interaction models are compared with the corresponding base model using:
- Likelihood-ratio test (LRT), with `p < 0.05`  
- Akaike information criterion (AIC), requiring an improvement of at least 2 points (`ΔAIC ≥ 2`)  

Interactions meeting both criteria are retained in the final model.

The script also:
- Reports regression coefficients, odds ratios, 95% confidence intervals, and p-values  
- Calculates 10-fold out-of-fold (OOF) AUC for the base and selected final models  
- Uses fold-specific scaling during cross-validation to avoid data leakage  
- Calculates model-based predicted probabilities for retained interactions  
- Evaluates continuous moderators at their 25th, 50th, and 75th percentiles  
- Evaluates predicted probabilities at −2 and +2 SD of the original PGS distribution  
- Estimates simple PGS slopes using `emmeans::emtrends()`  
- Tests pairwise contrasts between PGS slopes across moderator levels  

The simple-slope analyses assess whether the PGS-outcome association differs from zero at specific moderator levels, whereas slope contrasts assess whether these conditional PGS associations differ from one another.

---

### 5. Multiple testing correction

- `06_FDR_correction_multiple_testing.R`  

  Adjusts the relevant frequentist p-values using the Benjamini–Hochberg false discovery rate (FDR).

Multiple testing is controlled across the five outcome-specific tests within each PGS.

Statistical significance is defined as:

```text
q < 0.05 (two-sided)
```

Post-hoc simple-slope and slope-contrast analyses are used to interpret retained interaction effects and are not included in the primary FDR correction.

---

### 6. Bayesian sensitivity analyses

- `07_run_bayesian_sensitivity_analyses.R`  

  Re-estimates the models selected in the frequentist analysis using Bayesian logistic regression.

The Bayesian analysis does not repeat model selection. Instead, it uses the same:
- PGS-outcome combinations  
- Covariates  
- Complete-case samples  
- Site adjustments  
- Interaction terms retained by the frequentist model-selection procedure  

Regression coefficients are assigned weakly informative normal priors:

```text
Normal(0, 0.5)
```

on the log-odds scale, with autoscaling disabled. These priors regularise coefficient estimates toward zero when effects are weakly supported by the data.

Outputs include:
- Posterior coefficient estimates  
- Posterior odds ratios  
- 95% credible intervals  
- Posterior predicted probabilities for retained interaction models  

The Bayesian models are used as sensitivity analyses to assess whether the direction and magnitude of the selected frequentist associations remain stable after coefficient regularisation.

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

- GWAS summary statistics used to derive the latent-factor PGSs are publicly available from the corresponding source studies and should be downloaded separately.
- Individual-level cohort data are not included in this repository.

All scripts require user-defined paths and variable names to match the dataset being analysed.

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
- Subject identifiers should be harmonised across clinical, genetic, PGS, and covariate datasets before merging.  
- Continuous predictors are scaled before regression modelling.  
- Logistic regression models are fitted using complete-case data for the variables required by each PGS-outcome analysis.  
- The same complete-case sample is used when comparing the base and interaction models for a given PGS-outcome pair.  
- Site covariates can be specified separately for each outcome when required by the structure of the cohort.  
- Variable names, outcome definitions, moderator variables, and covariate structures should be adapted to the dataset being analysed.  
- Bayesian models reproduce the model specifications selected in the frequentist analysis and are intended as sensitivity analyses rather than a separate model-selection procedure.  

---

## Contact

For questions or collaboration, please contact the repository owner.
````
