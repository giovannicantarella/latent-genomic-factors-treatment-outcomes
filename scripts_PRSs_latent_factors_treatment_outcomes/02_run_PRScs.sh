#!/bin/bash
set -euo pipefail

# ============================================================
# Run PRS-CS across chromosomes using tmux
# ============================================================
# This script runs PRS-CS separately for chromosomes 1-22
# in parallel tmux sessions.
#
# Requirements:
# - Python 3 available
# - PRS-CS installed
# - tmux installed
# - LD reference panel available
#
# Notes:
# - One tmux session is created for each chromosome
# - Output logs are saved as log_chr<chromosome>.txt
# - The sample size should match the GWAS used to derive
#   the summary statistics
# ============================================================

# Insert your paths here
PRSCSPY="/path/to/PRScs.py"
REF_DIR="/path/to/ldblk_1kg_eur"
BIM_PREFIX="/path/to/target_dataset_prefix"
SST_FILE="/path/to/munged_or_filtered_summary_statistics.txt"
OUT_DIR="/path/to/output_directory"

OUT_PREFIX="${OUT_DIR}/output"

# Insert model parameters here
A=1.0
B=0.5
PHI=0.1

# Insert the appropriate GWAS sample size here
N_GWAS=71579

# MCMC settings
N_ITER=1000
N_BURNIN=500
THIN=5

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# Remove previous PRS-CS tmux sessions, if present
tmux ls -F '#S' 2>/dev/null | grep '^PRScs_chr' | xargs -r -n1 tmux kill-session -t || true

for chr in {1..22}; do
  SESSION_NAME="PRScs_chr${chr}"

  tmux new-session -d -s "${SESSION_NAME}" bash -lc "
    set -euo pipefail
    cd '${OUT_DIR}'

    PYTHONUNBUFFERED=1 python3 -u '${PRSCSPY}' \
      --ref_dir='${REF_DIR}' \
      --bim_prefix='${BIM_PREFIX}' \
      --sst_file='${SST_FILE}' \
      --a=${A} \
      --b=${B} \
      --phi=${PHI} \
      --n_gwas=${N_GWAS} \
      --n_iter=${N_ITER} \
      --n_burnin=${N_BURNIN} \
      --thin=${THIN} \
      --out_dir='${OUT_PREFIX}' \
      --chrom=\$(printf '%d' ${chr}) \
      --write_pst=FALSE \
      --write_psi=FALSE \
      2>&1 | tee log_chr${chr}.txt
  "
done

# Display active tmux sessions
tmux ls