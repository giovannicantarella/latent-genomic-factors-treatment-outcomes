#!/bin/bash
set -euo pipefail

# ============================================================
# Munge GWAS summary statistics using LDSC
# ============================================================
# This script prepares GWAS summary statistics for downstream
# analyses by:
# - filtering variants by MAF and INFO
# - retaining HapMap3 SNPs only
# - formatting effect sizes for signed statistics
#
# Requirements:
# - LDSC installed
# - Python 2 available
#
# Notes:
# - The resulting munged file can be used as input for PRS-CS
# - BETA is treated as the signed summary statistic
# ============================================================

# Insert your paths here
LDSC_DIR="/path/to/ldsc"
SUMSTATS_FILE="/path/to/input_summary_statistics.txt"
MERGE_ALLELES_FILE="/path/to/w_hm3.snplist"
OUT_PREFIX="/path/to/output_prefix"

# Insert the appropriate GWAS sample size here
N_GWAS=71579

python2 "${LDSC_DIR}/munge_sumstats.py" \
  --sumstats "${SUMSTATS_FILE}" \
  --out "${OUT_PREFIX}" \
  --N "${N_GWAS}" \
  --merge-alleles "${MERGE_ALLELES_FILE}" \
  --chunksize 500000 \
  --maf-min 0.01 \
  --info-min 0.9 \
  --signed-sumstats BETA,0 \
  --keep-maf