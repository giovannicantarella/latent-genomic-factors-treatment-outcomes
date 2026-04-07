# ============================================================
# Adjust main PRS effects for multiple testing using
# Benjamini-Hochberg FDR
# ============================================================
# Multiple testing is controlled across the four outcome-specific
# tests for each PRS variable.
#
# Significance threshold:
# - q < 0.10 (two-sided)
# ============================================================

results_df$p_fdr_bh <- NA_real_

for (prs_name in prs_vars) {
  idx <- which(results_df$prs == prs_name & results_df$term == prs_name)
  
  if (length(idx) > 0) {
    results_df$p_fdr_bh[idx] <- p.adjust(results_df$p[idx], method = "BH")
  }
}

results_df$significant_fdr_q10 <- ifelse(
  !is.na(results_df$p_fdr_bh) & results_df$p_fdr_bh < 0.10,
  1,
  0
)