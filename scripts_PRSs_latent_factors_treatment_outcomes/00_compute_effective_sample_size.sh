#!/bin/bash
set -euo pipefail

# ============================================================
# Compute mean effective sample size from GWAS summary statistics
# ============================================================
# This script calculates the mean effective sample size (Neff)
# from a tab-delimited GWAS summary statistics file.
#
# Formula:
#   Neff = 4 / (2 * MAF * (1 - MAF) * SE^2)
#
# Assumptions:
# - The input file is tab-delimited
# - The file contains a header row
# - MAF is in column 4
# - SE is in column 11
#
# Notes:
# - Quotation marks are removed before parsing values
# - Rows with invalid or missing MAF/SE values are excluded
#
# References:
# - https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics
# - https://isgw-forum.colorado.edu/t/backing-out-effective-sample-size/523
# ============================================================

# Usage:
#   bash 00_compute_effective_sample_size.sh /path/to/input_file.dat
#
# If no input file is provided, a placeholder filename is used.

FILE=${1:-"/path/to/input_file.dat"}

LC_NUMERIC=C awk -F'\t' '
NR > 1 {
  maf_raw = $4
  gsub(/"/, "", maf_raw)

  se_raw = $11
  gsub(/"/, "", se_raw)

  maf = maf_raw + 0
  se  = se_raw + 0

  if (maf <= 0 || maf >= 1) next
  if (se <= 0) next

  pq = 2 * maf * (1 - maf)
  if (pq <= 0) next

  neff = 4 / (pq * se * se)

  sum += neff
  count++
}
END {
  if (count > 0) {
    print "Valid rows:", count
    print "Mean Neff =", sum / count
  } else {
    print "No valid rows available to compute Neff."
  }
}
' "$FILE"