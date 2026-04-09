# ============================================================
# Project: Latent genomic factors and treatment outcomes
# Script: Logistic regression models for PRS association analyses
#
# Purpose:
# - read a prepared cohort dataset
# - fit logistic regression models across multiple outcomes and
#   polygenic score variables
# - test PRS main effects together with interactions involving
#   selected moderators
# - save a single log file for the full run
# - save a combined results table with beta, SE, z, p, odds ratios,
#   and 95% Wald confidence intervals
#
# Notes:
# - this script is written as a generic template and should be
#   adapted to the cohort being analysed
# - the input dataset is expected to be the output of a previous
#   preparation step
# - update file paths, outcome names, PRS variables, moderators,
#   and covariates before running
# - the same model structure is applied repeatedly across all
#   combinations of outcomes and PRS variables
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# 1. Input file and settings
# ============================================================
# Replace these placeholders with the input dataset and output
# directory for the cohort of interest.
#
# The prepared dataset is expected to include:
# - binary outcomes coded as 0/1
# - PRS variables already merged into the cohort dataset
# - optional moderator variables
# - site covariates and ancestry principal components
# ============================================================

input_file <- "/path/to/cohort_prepared_with_PRS.csv"
output_dir <- "/path/to/output_directory"

# Outcomes to be analysed.
# These should be binary variables already present in the dataset.
outcomes <- c(
  "OUTCOME_1",
  "OUTCOME_2",
  "OUTCOME_3"
)

# PRS variables to be tested.
# These are assumed to be standardised or otherwise comparable
# across participants before entering the models.
prs_vars <- c(
  "PRS_A_z",
  "PRS_B_z",
  "PRS_C_z",
  "PRS_D_z",
  "PRS_E_z"
)

# Site and ancestry covariates.
# Adjust these vectors to reflect the covariates available in the
# cohort being analysed.
cov_sites <- paste0("site", 1:9)
cov_pcs   <- paste0("PC", 1:9)

# Moderators to be tested together with each PRS.
# The first moderator below is treated as a continuous clinical
# severity variable in the plotting functions.
severity_var <- "BASELINE_SEVERITY"
moderators <- c(severity_var, "age", "sex")

# Threshold used to decide whether a plot should be generated for
# a given outcome–PRS combination.
p_plot_threshold <- 0.10

# ============================================================
# 2. Helper functions
# ============================================================
# The functions below are used repeatedly across the script to:
# - check outcome coding
# - extract model summaries in a compact format
# - build reference newdata objects for plotting
# - save plots for main effects and interactions
# ============================================================

# Convert an outcome to 0/1 and stop if values outside 0/1 are found.
# This makes the expected input explicit before entering the model loop.
as_binary_01 <- function(x, name = "outcome") {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) x <- suppressWarnings(as.numeric(x))
  
  ux <- sort(unique(x[!is.na(x)]))
  if (!all(ux %in% c(0, 1))) {
    stop(
      name, " is not coded as 0/1. Observed values: ",
      paste(ux, collapse = ", ")
    )
  }
  
  as.integer(x)
}

# Extract coefficient estimates from glm() output and add
# Wald-based odds ratios and confidence intervals.
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

# Print a compact version of a coefficient table.
# This is useful for the log file, where keeping the output
# readable matters more than preserving many decimals.
print_table <- function(df, title = NULL) {
  if (!is.null(title)) cat("\n", title, "\n", sep = "")
  
  df_print <- df
  num_cols <- c("estimate", "se", "z", "p", "OR", "CI_low", "CI_high")
  
  for (cc in intersect(num_cols, names(df_print))) {
    df_print[[cc]] <- signif(df_print[[cc]], 4)
  }
  
  print(df_print, row.names = FALSE)
}

# Create a baseline newdata object used for plotting model-based
# predicted probabilities. Continuous covariates are fixed at their
# mean; binary or factor-like terms are fixed at a default or common
# reference value.
make_baseline_newdata <- function(dat, prs, cov_sites, cov_pcs,
                                  severity_var,
                                  moderator = NULL,
                                  grid_vals = NULL) {
  nd <- list()
  
  for (s in cov_sites) {
    if (s %in% names(dat)) nd[[s]] <- 0
  }
  
  for (pc in cov_pcs) {
    if (pc %in% names(dat)) nd[[pc]] <- mean(dat[[pc]], na.rm = TRUE)
  }
  
  if ("age" %in% names(dat)) {
    nd[["age"]] <- mean(dat[["age"]], na.rm = TRUE)
  }
  
  if (severity_var %in% names(dat)) {
    nd[[severity_var]] <- mean(dat[[severity_var]], na.rm = TRUE)
  }
  
  if ("sex" %in% names(dat)) {
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
  }
  
  nd[[prs]] <- 0
  
  if (!is.null(moderator) && !is.null(grid_vals)) {
    nd[[moderator]] <- grid_vals
  }
  
  as.data.frame(nd, stringsAsFactors = FALSE)
}

# Plot predicted probability across the observed PRS range while
# holding the remaining covariates fixed.
save_plot_main_prs <- function(fit, dat, outcome, prs, outdir,
                               cov_sites, cov_pcs, severity_var) {
  rng <- quantile(dat[[prs]], probs = c(0.05, 0.95), na.rm = TRUE)
  xseq <- seq(rng[1], rng[2], length.out = 120)
  
  nd <- make_baseline_newdata(
    dat = dat,
    prs = prs,
    cov_sites = cov_sites,
    cov_pcs = cov_pcs,
    severity_var = severity_var
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

# Plot a PRS × continuous moderator interaction.
# Predicted probability is shown across the moderator range for
# three representative PRS values (-1, 0, +1).
save_plot_interaction_cont <- function(fit, dat, outcome, prs, moderator,
                                       outdir, cov_sites, cov_pcs, severity_var) {
  rng <- quantile(dat[[moderator]], probs = c(0.05, 0.95), na.rm = TRUE)
  mseq <- seq(rng[1], rng[2], length.out = 120)
  
  prs_levels <- c(-1, 0, 1)
  file_out <- file.path(
    outdir,
    paste0(outcome, "__", prs, "__x__", gsub("\\.", "_", moderator), ".png")
  )
  
  png(file_out, width = 1400, height = 900, res = 180)
  
  plot(
    NA,
    xlim = range(mseq),
    ylim = c(0, 1),
    xlab = moderator,
    ylab = paste0("P(", outcome, " = 1)"),
    main = paste0(outcome, " | ", prs, " × ", moderator, " (PRS = -1 / 0 / +1)")
  )
  
  for (lv in prs_levels) {
    nd <- make_baseline_newdata(
      dat = dat,
      prs = prs,
      cov_sites = cov_sites,
      cov_pcs = cov_pcs,
      severity_var = severity_var,
      moderator = moderator,
      grid_vals = mseq
    )
    nd <- nd[rep(1, length(mseq)), , drop = FALSE]
    nd[[moderator]] <- mseq
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

# Plot a PRS × binary moderator interaction.
save_plot_interaction_binary <- function(fit, dat, outcome, prs, moderator,
                                         outdir, cov_sites, cov_pcs, severity_var) {
  prs_levels <- c(-1, 0, 1)
  moderator_levels <- c(0, 1)
  
  file_out <- file.path(
    outdir,
    paste0(outcome, "__", prs, "__x__", moderator, ".png")
  )
  png(file_out, width = 1200, height = 900, res = 180)
  
  plot(
    NA,
    xlim = c(0.5, 3.5),
    ylim = c(0, 1),
    xaxt = "n",
    xlab = "PRS level",
    ylab = paste0("P(", outcome, " = 1)"),
    main = paste0(outcome, " | ", prs, " × ", moderator)
  )
  axis(1, at = 1:3, labels = paste0(prs_levels))
  
  x <- 1:3
  for (lv in moderator_levels) {
    nd <- make_baseline_newdata(
      dat = dat,
      prs = prs,
      cov_sites = cov_sites,
      cov_pcs = cov_pcs,
      severity_var = severity_var
    )
    nd <- nd[rep(1, length(prs_levels)), , drop = FALSE]
    nd[[prs]] <- prs_levels
    nd[[moderator]] <- lv
    
    pr <- predict(fit, newdata = nd, type = "response")
    lines(x, pr, type = "b", lwd = 2, pch = 16)
  }
  
  legend(
    "topright",
    legend = paste0(moderator, " = ", moderator_levels),
    lwd = 2,
    pch = 16,
    bty = "n"
  )
  
  dev.off()
  
  cat("Saved plot: ", file_out, "\n", sep = "")
}

# Extract the p-value for a specific model term if present.
get_term_pvalue <- function(tab, term_name) {
  idx <- which(tab$term == term_name)
  if (length(idx) == 0) return(NA_real_)
  tab$p[idx[1]]
}

# ============================================================
# 3. Read input data
# ============================================================
# The input dataset should already contain all variables required
# for modelling, including outcomes, PRSs, moderators, sites,
# and principal components.
# ============================================================

cohort_df <- read.csv(input_file, stringsAsFactors = FALSE)

# ============================================================
# 4. Safety checks
# ============================================================
# Before fitting models, check that all common covariates,
# outcomes, moderators, and PRS variables are present.
# These checks fail early and make debugging easier.
# ============================================================

needed_common <- unique(c(cov_sites, cov_pcs, moderators))
missing_common <- setdiff(needed_common, names(cohort_df))
if (length(missing_common) > 0) {
  stop(
    "Missing common covariates in the cohort dataset: ",
    paste(missing_common, collapse = ", ")
  )
}

missing_outcomes <- setdiff(outcomes, names(cohort_df))
if (length(missing_outcomes) > 0) {
  stop(
    "Missing outcomes in the cohort dataset: ",
    paste(missing_outcomes, collapse = ", ")
  )
}

missing_prs <- setdiff(prs_vars, names(cohort_df))
if (length(missing_prs) > 0) {
  stop(
    "Missing PRS variables in the cohort dataset: ",
    paste(missing_prs, collapse = ", ")
  )
}

# ============================================================
# 5. Output directories and log files
# ============================================================
# Each run is saved in a timestamped folder. This helps keep
# repeated analyses separate and avoids overwriting previous
# logs, results, and figures.
# ============================================================

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(output_dir, paste0("cohort_glm_prs_", timestamp))
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- file.path(run_dir, "plots_p_lt_threshold")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(run_dir, paste0("cohort_glm_prs_", timestamp, ".txt"))
csv_file <- file.path(run_dir, paste0("cohort_glm_prs_results_", timestamp, ".csv"))

# ============================================================
# 6. Main analysis loop
# ============================================================
# The same model structure is fitted for each combination of:
# - outcome
# - PRS variable
#
# The regression includes:
# - PRS main effect
# - PRS × severity interaction
# - PRS × age interaction
# - PRS × sex interaction
# - site covariates
# - ancestry principal components
#
# All console output is redirected to a single log file while
# still being printed on screen.
# ============================================================

all_results <- list()
result_index <- 0

sink(log_file, split = TRUE)

cat("=== Logistic regression loop for PRS analyses started ===\n")
cat("Timestamp: ", timestamp, "\n", sep = "")
cat("Run directory: ", run_dir, "\n", sep = "")
cat("Plot threshold: p < ", p_plot_threshold, "\n", sep = "")

for (outcome in outcomes) {
  
  cat("\n\n============================================================\n")
  cat("OUTCOME: ", outcome, "\n", sep = "")
  cat("Raw class balance:\n")
  print(table(cohort_df[[outcome]], useNA = "ifany"))
  
  # Recode outcome to 0/1 explicitly before model fitting.
  cohort_df[[outcome]] <- as_binary_01(cohort_df[[outcome]], name = outcome)
  
  for (prs in prs_vars) {
    
    cat("\n------------------------------------------------------------\n")
    cat("PRS: ", prs, "\n", sep = "")
    
    # Construct the right-hand side of the model formula.
    # This uses the same covariate structure for each PRS.
    rhs <- paste(
      paste0(prs, "*", severity_var),
      paste0(prs, "*age"),
      paste0(prs, "*sex"),
      paste(c(cov_sites, cov_pcs), collapse = " + "),
      sep = " + "
    )
    
    formula_obj <- as.formula(paste(outcome, "~", rhs))
    
    fit <- glm(formula_obj, data = cohort_df, family = binomial())
    
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
    
    # Identify the PRS-related terms that control the plotting logic.
    term_main      <- prs
    term_severity  <- paste0(prs, ":", severity_var)
    term_age       <- paste0(prs, ":age")
    term_sex       <- paste0(prs, ":sex")
    
    p_main     <- get_term_pvalue(tab, term_main)
    p_severity <- get_term_pvalue(tab, term_severity)
    p_age      <- get_term_pvalue(tab, term_age)
    p_sex      <- get_term_pvalue(tab, term_sex)
    
    do_plot <- any(c(p_main, p_severity, p_age, p_sex) < p_plot_threshold, na.rm = TRUE)
    
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
          cov_pcs = cov_pcs,
          severity_var = severity_var
        )
      }
      
      if (!is.na(p_severity) && p_severity < p_plot_threshold) {
        save_plot_interaction_cont(
          fit = fit,
          dat = dat_used,
          outcome = outcome,
          prs = prs,
          moderator = severity_var,
          outdir = plot_dir,
          cov_sites = cov_sites,
          cov_pcs = cov_pcs,
          severity_var = severity_var
        )
      }
      
      if (!is.na(p_age) && p_age < p_plot_threshold) {
        save_plot_interaction_cont(
          fit = fit,
          dat = dat_used,
          outcome = outcome,
          prs = prs,
          moderator = "age",
          outdir = plot_dir,
          cov_sites = cov_sites,
          cov_pcs = cov_pcs,
          severity_var = severity_var
        )
      }
      
      if (!is.na(p_sex) && p_sex < p_plot_threshold) {
        save_plot_interaction_binary(
          fit = fit,
          dat = dat_used,
          outcome = outcome,
          prs = prs,
          moderator = "sex",
          outdir = plot_dir,
          cov_sites = cov_sites,
          cov_pcs = cov_pcs,
          severity_var = severity_var
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
# All coefficient tables are bound together in a single dataset
# and saved as one csv file for downstream inspection or reporting.
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