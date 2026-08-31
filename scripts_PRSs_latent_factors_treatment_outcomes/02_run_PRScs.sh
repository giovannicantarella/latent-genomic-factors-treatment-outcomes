#!/bin/bash
set -euo pipefail
# ============================================================
# Run PRS-CS across chromosomes using tmux
# ============================================================
# Script: Run PRS-CS across chromosomes using tmux
#
# Purpose:
# - run PRS-CS separately for chromosomes 1-22
# - launch one parallel tmux session per chromosome
# - save one log file for each chromosome-specific run
#
# Workflow:
# - PRS-CS is executed independently for each chromosome
# - each chromosome is sent to a dedicated tmux session
# - logs are written to chromosome-specific text files
# - output files can later be combined downstream if needed
#
# Requirements:
# - Python 3 available in the environment
# - PRS-CS installed and callable via the selected Python interpreter
# - tmux installed
# - LD reference panel available
# - target dataset .bim/.bed/.fam prefix available
#
# Notes:
# - this script is written as a generic template and should be
#   adapted to the dataset and GWAS being used
# - one tmux session is created for each chromosome
# - output logs are saved as log_chr<chromosome>.txt
# - the GWAS sample size should match the summary statistics
#   used as input to PRS-CS
# - the output prefix is shared across chromosomes; PRS-CS will
#   append chromosome-specific information to the generated files
# ============================================================

# ============================================================
# 1. Input paths
# ============================================================
# Replace these placeholders with the relevant paths for the
# current analysis.
#
# PRSCSPY   : path to PRScs.py
# REF_DIR   : directory containing the LD reference panel
# BIM_PREFIX: prefix of the target genotype dataset
#            (without .bed/.bim/.fam extension)
# SST_FILE  : munged or filtered summary statistics
# OUT_DIR   : directory where logs and PRS-CS outputs will be saved
# ============================================================

PRSCSPY="/path/to/PRScs.py"
REF_DIR="/path/to/ldblk_1kg_eur"
BIM_PREFIX="/path/to/target_dataset_prefix"
SST_FILE="/path/to/munged_or_filtered_summary_statistics.txt"
OUT_DIR="/path/to/output_directory"

# Output prefix used by PRS-CS.
# This is not a directory name but the prefix used to construct
# the output filenames written by PRS-CS.
OUT_PREFIX="${OUT_DIR}/output"

# ============================================================
# 2. Model parameters
# ============================================================
# These are the main PRS-CS hyperparameters.
# Update them according to the analytical plan being used.
#
# A, B, PHI:
# - control the prior specification in PRS-CS
# - PHI is often the parameter most frequently tuned by the user
# ============================================================

A=1.0
B=0.5
PHI=0.1

# GWAS sample size used to derive the summary statistics.
# This should correspond to the discovery GWAS underlying SST_FILE.
N_GWAS=71579

# ============================================================
# 3. MCMC settings
# ============================================================
# These settings control the PRS-CS MCMC procedure.
#
# N_ITER  : total number of iterations
# N_BURNIN: burn-in iterations discarded before posterior summaries
# THIN    : thinning interval
#
# Increase these values if a more stable posterior estimate is needed,
# keeping in mind the additional computational cost.
# ============================================================

N_ITER=1000
N_BURNIN=500
THIN=5

# ============================================================
# 4. Prepare output directory
# ============================================================
# The output directory is created if it does not already exist.
# The script then moves into that directory so that log files are
# written locally and are easier to inspect.
# ============================================================

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# ============================================================
# 5. Remove previous PRS-CS tmux sessions
# ============================================================
# If old chromosome-specific tmux sessions from a previous run are
# still active, they are removed here to avoid conflicts with the
# new run. Sessions not matching the chosen naming convention are
# left untouched.
# ============================================================

tmux ls -F '#S' 2>/dev/null | grep '^PRScs_chr' | xargs -r -n1 tmux kill-session -t || true

# ============================================================
# 6. Launch one tmux session per chromosome
# ============================================================
# PRS-CS is run separately for chromosomes 1 through 22.
# Each run is started in detached mode so that all chromosome jobs
# can proceed in parallel in the background.
#
# For each chromosome:
# - a session named PRScs_chr<chr> is created
# - PRScs.py is called with chromosome-specific arguments
# - stdout and stderr are both saved to log_chr<chr>.txt
# ============================================================

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

# ============================================================
# 7. Display active tmux sessions
# ============================================================
# At the end of the script, list active tmux sessions so the user
# can quickly verify that chromosome-specific jobs have started.
#
# Useful commands after launch:
# - tmux attach -t PRScs_chr1
# - tmux capture-pane -pt PRScs_chr1
# - tmux kill-session -t PRScs_chr1
# ============================================================

tmux ls
