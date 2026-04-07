#!/bin/bash
set -euo pipefail

# ============================================================
# Combine PRS-CS weights and compute polygenic scores
# ============================================================
# This script:
# - merges chromosome-specific PRS-CS weight files
# - computes polygenic scores in the target dataset using PLINK2
#
# Requirements:
# - PLINK2 installed
#
# Notes:
# - The chromosome-specific files are assumed to follow the
#   standard PRS-CS naming convention
# - The column indices in --score should match the format
#   of the PRS-CS output files
# ============================================================

# Insert your paths here
WEIGHTS_PREFIX="/path/to/output_pst_eff_a1_b0.5_phi1e-01_chr"
COMBINED_WEIGHTS_FILE="/path/to/PRScs_weights_all.txt"
BFILE_PREFIX="/path/to/target_dataset_prefix"
OUT_PREFIX="/path/to/PRScs_scores"

# Combine chromosome-specific weight files
head -n 1 "${WEIGHTS_PREFIX}1.txt" > "${COMBINED_WEIGHTS_FILE}"
tail -n +2 -q "${WEIGHTS_PREFIX}"{1..22}.txt >> "${COMBINED_WEIGHTS_FILE}"

# Compute polygenic scores
plink2 \
  --bfile "${BFILE_PREFIX}" \
  --score "${COMBINED_WEIGHTS_FILE}" 2 4 6 no-mean-imputation \
  --out "${OUT_PREFIX}" \
  --threads 4 \
  --memory 24000