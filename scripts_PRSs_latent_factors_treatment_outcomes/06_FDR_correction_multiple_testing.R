# ============================================================
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Adjust selected PGS effects for multiple testing using
# Benjamini-Hochberg FDR
# ============================================================
# Multiple testing is controlled across the five outcome-specific
# tests for each PGS variable.
#
# For each PGS-outcome pair, the focal effect corresponds to:
# - the PGS main effect for models without interactions;
# - the PGS-by-moderator term for models with an interaction.
#
# Significance threshold:
# - q < 0.05 (two-sided)
# ============================================================

results_df$p_fdr_bh <- NA_real_

for (pgs_name in pgs_vars) {
  
  pgs_rows <- which(results_df$pgs == pgs_name)
  
  focal_idx <- unlist(
    lapply(
      split(pgs_rows, results_df$outcome[pgs_rows]),
      function(idx) {
        
        interaction_idx <- idx[
          grepl(":", results_df$term[idx], fixed = TRUE)
        ]
        
        main_idx <- idx[
          results_df$term[idx] == pgs_name
        ]
        
        if (length(interaction_idx) == 1) {
          interaction_idx
        } else if (length(interaction_idx) == 0 &&
                   length(main_idx) == 1) {
          main_idx
        } else {
          stop(
            paste(
              "A unique focal effect could not be identified for",
              pgs_name,
              "and outcome",
              unique(results_df$outcome[idx])
            )
          )
        }
      }
    )
  )
  
  if (length(focal_idx) != 5) {
    stop(
      paste(
        "Expected five outcome-specific tests for",
        pgs_name,
        "but found",
        length(focal_idx)
      )
    )
  }
  
  results_df$p_fdr_bh[focal_idx] <- p.adjust(
    results_df$p[focal_idx],
    method = "BH"
  )
}

results_df$significant_fdr_q05 <- ifelse(
  !is.na(results_df$p_fdr_bh) &
    results_df$p_fdr_bh < 0.05,
  1,
  0
)
