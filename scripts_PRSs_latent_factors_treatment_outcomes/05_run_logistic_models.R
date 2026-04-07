# ============================================================
# Run logistic regression models for GSRD PRS analyses
# ============================================================
# This script:
# - reads a prepared GSRD dataset
# - fits logistic regression models across multiple outcomes
#   and PRS variables
# - includes PRS interactions with selected moderators
# - writes a single log file
# - saves a combined results table with beta, SE, z, p,
#   odds ratios, and 95% Wald confidence intervals
# - generates plots only when at least one PRS-related term
#   reaches the specified p-value threshold
#
# Notes:
# - Update file paths before running
# - The input dataset is expected to be the output of the
#   preparation script
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# 1. Input file and settings
# ============================================================

# Insert your paths here
input_file  <- "/path/to/GSRD_prepared_with_PRS.csv"
output_dir  <- "/path/to/output_directory"

outcomes <- c(
  "remission",
  "response",
  "resistance",
  "augmentation_antipsych",
  "augmentation_moodstab"
)

prs_vars <- c(
  "PRS_ThoughtDisorders_z",
  "PRS_Compulsive_z",
  "PRS_Neuro_z",
  "PRS_SUD_z",
  "PRS_Internalizing_z"
)

cov_sites <- paste0("site", 1:9)
cov_pcs   <- paste0("PC", 1:9)
moderators <- c("MADRS.Retrospective", "age", "sex")

p_plot_threshold <- 0.10

# ============================================================
# 2. Helper functions
# ============================================================

as_binary_01 <- function(x, name = "outcome") {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) x <- suppressWarnings(as.numeric(x))
  
  ux <- sort(unique(x[!is.na(x)]))
  if (!all(ux %in% c(0, 1))) {
    stop(name, " is not coded as 0/1. Observed values: ",
         paste(ux, collapse = ", "))
  }
  
  as.integer(x)
}

tidy_glm_or <- function(fit) {
  sm <- summary(fit)$coefficients
  
  out <- data.frame(
    term = rownames(sm),
    estimate = sm[, "Estimate"],
    se = sm[, "Std. Error"],
    z = sm[, "z value"],
    p = sm[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  zcrit <- qnorm(0.975)
  out$ci_low_beta  <- out$estimate - zcrit * out$se
  out$ci_high_beta <- out$estimate + zcrit * out$se
  out$OR      <- exp(out$estimate)
  out$CI_low  <- exp(out$ci_low_beta)
  out$CI_high <- exp(out$ci_high_beta)
  
  out
}

print_table <- function(df, title = NULL) {
  if (!is.null(title)) cat("\n", title, "\n", sep = "")
  
  df_print <- df
  num_cols <- c("estimate", "se", "z", "p", "OR", "CI_low", "CI_high")
  
  for (cc in intersect(num_cols, names(df_print))) {
    df_print[[cc]] <- signif(df_print[[cc]], 4)
  }
  
  print(df_print, row.names = FALSE)
}

make_baseline_newdata <- function(dat, prs, cov_sites, cov_pcs, mod = NULL, grid_vals = NULL) {
  nd <- list()
  
  for (s in cov_sites) {
    nd[[s]] <- 0
  }
  
  for (pc in cov_pcs) {
    nd[[pc]] <- mean(dat[[pc]], na.rm = TRUE)
  }
  
  nd[["age"]] <- mean(dat[["age"]], na.rm = TRUE)
  nd[["MADRS.Retrospective"]] <- mean(dat[["MADRS.Retrospective"]], na.rm = TRUE)
  
  sx <- dat[["sex"]]
  if (is.factor(sx)) sx <- as.character(sx)
  if (is.character(sx)) sx <- suppressWarnings(as.numeric(sx))
  ux <- sort(unique(sx[!is.na(sx)]))
  
  if (all(ux %in% c(0, 1))) {
    nd[["sex"]] <- 0
  } else {
    tab <- sort(table(dat[["sex"]]), decreasing = TRUE)
    nd[["sex"]] <- names(tab)[1]
  }
  
  nd[[prs]] <- 0
  
  if (!is.null(mod) && !is.null(grid_vals)) {
    nd[[mod]] <- grid_vals
  }
  
  as.data.frame(nd, stringsAsFactors = FALSE)
}

save_plot_main_prs <- function(fit, dat, outcome, prs, outdir, cov_sites, cov_pcs) {
  rng <- quantile(dat[[prs]], probs = c(0.05, 0.95), na.rm = TRUE)
  xseq <- seq(rng[1], rng[2], length.out = 120)
  
  nd <- make_baseline_newdata(
    dat = dat,
    prs = prs,
    cov_sites = cov_sites,
    cov_pcs = cov_pcs
  )
  nd <- nd[rep(1, length(xseq)), , drop = FALSE]
  nd[[prs]] <- xseq
  
  pr <- predict(fit, newdata = nd, type = "response")
  
  file_out <- file.path(outdir, paste0(outcome, "__", prs, "__mainPRS.png"))
  png(file_out, width = 1400, height = 900, res = 180)
  plot(
    xseq, pr,
    type = "l",
    xlab = prs,
    ylab = paste0("P(", outcome, " = 1)"),
    main = paste0(outcome, " | ", prs, " (main effect)")
  )
  abline(h = 0.5, lty = 3)
  dev.off()
  
  cat("Saved plot: ", file_out, "\n", sep = "")
}

save_plot_interaction_cont <- function(fit, dat, outcome, prs, mod, outdir, cov_sites, cov_pcs) {
  rng <- quantile(dat[[mod]], probs = c(0.05, 0.95), na.rm = TRUE)
  mseq <- seq(rng[1], rng[2], length.out = 120)
  
  prs_levels <- c(-1, 0, 1)
  file_out <- file.path(
    outdir,
    paste0(outcome, "__", prs, "__x__", gsub("\\.", "_", mod), ".png")
  )
  
  png(file_out, width = 1400, height = 900, res = 180)
  
  plot(
    NA,
    xlim = range(mseq),
    ylim = c(0, 1),
    xlab = mod,
    ylab = paste0("P(", outcome, " = 1)"),
    main = paste0(outcome, " | ", prs, " × ", mod, " (PRS = -1 / 0 / +1)")
  )
  
  for (lv in prs_levels) {
    nd <- make_baseline_newdata(
      dat = dat,
      prs = prs,
      cov_sites = cov_sites,
      cov_pcs = cov_pcs,
      mod = mod,
      grid_vals = mseq
    )
    nd <- nd[rep(1, length(mseq)), , drop = FALSE]
    nd[[mod]] <- mseq
    nd[[prs]] <- lv
    
    pr <- predict(fit, newdata = nd, type = "response")
    lines(mseq, pr, lwd = 2)
  }
  
  legend(
    "topright",
    legend = paste0(prs, " = ", prs_levels),
    lwd = 2,
    bty = "n"
  )
  
  dev.off()
  
  cat("Saved plot: ", file_out, "\n", sep = "")
}

save_plot_interaction_sex <- function(fit, dat, outcome, prs, outdir, cov_sites, cov_pcs) {
  prs_levels <- c(-1, 0, 1)
  sexes <- c(0, 1)
  
  file_out <- file.path(outdir, paste0(outcome, "__", prs, "__x__sex.png"))
  png(file_out, width = 1200, height = 900, res = 180)
  
  plot(
    NA,
    xlim = c(0.5, 3.5),
    ylim = c(0, 1),
    xaxt = "n",
    xlab = "PRS level",
    ylab = paste0("P(", outcome, " = 1)"),
    main = paste0(outcome, " | ", prs, " × sex")
  )
  axis(1, at = 1:3, labels = paste0(prs_levels))
  
  x <- 1:3
  for (sx in sexes) {
    nd <- make_baseline_newdata(
      dat = dat,
      prs = prs,
      cov_sites = cov_sites,
      cov_pcs = cov_pcs
    )
    nd <- nd[rep(1, length(prs_levels)), , drop = FALSE]
    nd[[prs]] <- prs_levels
    nd[["sex"]] <- sx
    
    pr <- predict(fit, newdata = nd, type = "response")
    lines(x, pr, type = "b", lwd = 2, pch = 16)
  }
  
  legend(
    "topright",
    legend = c("sex = 0", "sex = 1"),
    lwd = 2,
    pch = 16,
    bty = "n"
  )
  
  dev.off()
  
  cat("Saved plot: ", file_out, "\n", sep = "")
}

get_term_pvalue <- function(tab, term_name) {
  idx <- which(tab$term == term_name)
  if (length(idx) == 0) return(NA_real_)
  tab$p[idx[1]]
}

# ============================================================
# 3. Read input data
# ============================================================

GSRD <- read.csv(input_file, stringsAsFactors = FALSE)

# ============================================================
# 4. Safety checks
# ============================================================

needed_common <- unique(c(cov_sites, cov_pcs, moderators))
missing_common <- setdiff(needed_common, names(GSRD))
if (length(missing_common) > 0) {
  stop("Missing common covariates in GSRD: ",
       paste(missing_common, collapse = ", "))
}

missing_outcomes <- setdiff(outcomes, names(GSRD))
if (length(missing_outcomes) > 0) {
  stop("Missing outcomes in GSRD: ",
       paste(missing_outcomes, collapse = ", "))
}

missing_prs <- setdiff(prs_vars, names(GSRD))
if (length(missing_prs) > 0) {
  stop("Missing PRS variables in GSRD: ",
       paste(missing_prs, collapse = ", "))
}

# ============================================================
# 5. Output directories and log files
# ============================================================

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(output_dir, paste0("GSRD_GLM_TDPRS_", timestamp))
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- file.path(run_dir, "plots_p_lt_0p10")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(run_dir, paste0("GSRD_GLM_TDPRS_", timestamp, ".txt"))
csv_file <- file.path(run_dir, paste0("GSRD_GLM_TDPRS_results_", timestamp, ".csv"))

# ============================================================
# 6. Main analysis loop
# ============================================================

all_results <- list()
result_index <- 0

sink(log_file, split = TRUE)

cat("=== GSRD logistic regression loop started ===\n")
cat("Timestamp: ", timestamp, "\n", sep = "")
cat("Run directory: ", run_dir, "\n", sep = "")
cat("Plot threshold: p < ", p_plot_threshold, "\n", sep = "")

for (outcome in outcomes) {
  
  cat("\n\n============================================================\n")
  cat("OUTCOME: ", outcome, "\n", sep = "")
  cat("Raw class balance:\n")
  print(table(GSRD[[outcome]], useNA = "ifany"))
  
  GSRD[[outcome]] <- as_binary_01(GSRD[[outcome]], name = outcome)
  
  for (prs in prs_vars) {
    
    cat("\n------------------------------------------------------------\n")
    cat("PRS: ", prs, "\n", sep = "")
    
    rhs <- paste(
      paste0(prs, "*MADRS.Retrospective"),
      paste0(prs, "*age"),
      paste0(prs, "*sex"),
      paste(c(cov_sites, cov_pcs), collapse = " + "),
      sep = " + "
    )
    
    formula_obj <- as.formula(paste(outcome, "~", rhs))
    
    fit <- glm(formula_obj, data = GSRD, family = binomial())
    
    n_used <- nobs(fit)
    y_used <- model.response(model.frame(fit))
    events <- sum(y_used == 1, na.rm = TRUE)
    
    cat("N used (complete cases): ", n_used, "\n", sep = "")
    cat("Events (", outcome, " = 1): ", events, "\n", sep = "")
    
    cat("\nsummary(glm):\n")
    print(summary(fit))
    
    tab <- tidy_glm_or(fit)
    tab$outcome <- outcome
    tab$prs <- prs
    tab$n_used <- n_used
    tab$events <- events
    
    result_index <- result_index + 1
    all_results[[result_index]] <- tab
    
    print_table(
      tab[, c(
        "outcome", "prs", "term", "estimate", "se", "z", "p",
        "OR", "CI_low", "CI_high", "n_used", "events"
      )],
      title = "Full coefficient table (Wald OR and 95% CI):"
    )
    
    term_main  <- prs
    term_madrs <- paste0(prs, ":MADRS.Retrospective")
    term_age   <- paste0(prs, ":age")
    term_sex   <- paste0(prs, ":sex")
    
    p_main  <- get_term_pvalue(tab, term_main)
    p_madrs <- get_term_pvalue(tab, term_madrs)
    p_age   <- get_term_pvalue(tab, term_age)
    p_sex   <- get_term_pvalue(tab, term_sex)
    
    do_plot <- any(c(p_main, p_madrs, p_age, p_sex) < p_plot_threshold, na.rm = TRUE)
    
    if (do_plot) {
      cat("\nPlot trigger: at least one PRS-related term has p < ",
          p_plot_threshold, "\n", sep = "")
      
      dat_used <- model.frame(fit)
      
      if (!is.na(p_main) && p_main < p_plot_threshold) {
        save_plot_main_prs(
          fit = fit,
          dat = dat_used,
          outcome = outcome,
          prs = prs,
          outdir = plot_dir,
          cov_sites = cov_sites,
          cov_pcs = cov_pcs
        )
      }
      
      if (!is.na(p_madrs) && p_madrs < p_plot_threshold) {
        save_plot_interaction_cont(
          fit = fit,
          dat = dat_used,
          outcome = outcome,
          prs = prs,
          mod = "MADRS.Retrospective",
          outdir = plot_dir,
          cov_sites = cov_sites,
          cov_pcs = cov_pcs
        )
      }
      
      if (!is.na(p_age) && p_age < p_plot_threshold) {
        save_plot_interaction_cont(
          fit = fit,
          dat = dat_used,
          outcome = outcome,
          prs = prs,
          mod = "age",
          outdir = plot_dir,
          cov_sites = cov_sites,
          cov_pcs = cov_pcs
        )
      }
      
      if (!is.na(p_sex) && p_sex < p_plot_threshold) {
        save_plot_interaction_sex(
          fit = fit,
          dat = dat_used,
          outcome = outcome,
          prs = prs,
          outdir = plot_dir,
          cov_sites = cov_sites,
          cov_pcs = cov_pcs
        )
      }
      
    } else {
      cat("\nNo plots generated: all PRS-related terms have p >= ",
          p_plot_threshold, " or are missing.\n", sep = "")
    }
  }
}

# ============================================================
# 7. Save combined results
# ============================================================

results_df <- do.call(rbind, all_results)

results_df <- results_df[, c(
  "outcome", "prs", "term", "estimate", "se", "z", "p",
  "OR", "CI_low", "CI_high", "n_used", "events",
  "ci_low_beta", "ci_high_beta"
)]

write.csv(results_df, csv_file, row.names = FALSE)

cat("\n\nSaved results table: ", csv_file, "\n", sep = "")
cat("Saved log file: ", log_file, "\n", sep = "")
cat("Plots directory: ", plot_dir, "\n", sep = "")
cat("\n=== DONE ===\n")

sink()