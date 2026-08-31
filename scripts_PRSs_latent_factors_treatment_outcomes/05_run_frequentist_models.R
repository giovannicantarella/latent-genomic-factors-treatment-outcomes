# ============================================================
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Frequentist logistic regression models for PGS association analyses
#
# Purpose:
# - fit logistic regression models for multiple binary outcomes
# - test associations between polygenic scores (PGSs) and outcomes
# - evaluate selected PGS interaction terms
# - compare nested models using likelihood-ratio tests and AIC
# - estimate out-of-fold AUC with fold-specific scaling
# - calculate predicted probabilities for retained interactions
# - estimate simple PGS slopes and contrasts between slopes
# - save analysis-ready summary tables
#
# Notes:
# - this script is written as a generic template and should be
#   adapted to the variable names and structure of the cohort
# - continuous variables are scaled using Gelman scaling:
#   (x - mean) / (2 * SD)
# - interaction models are retained when both:
#     LRT p < 0.05
#     delta AIC >= 2 relative to the base model
# - post-hoc interaction analyses are descriptive and intended
#   to aid interpretation of retained interaction effects
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(emmeans)
library(readr)
library(writexl)


# ============================================================
# 1. User-defined settings
# ============================================================

data_file <- "analysis_ready_dataset.csv"
output_dir <- "frequentist_outputs"

# Replace these variable names with those used in the cohort.
prs_vars <- c(
  "PRS_1",
  "PRS_2",
  "PRS_3",
  "PRS_4",
  "PRS_5"
)

outcomes <- c(
  "OUTCOME_1",
  "OUTCOME_2",
  "OUTCOME_3",
  "OUTCOME_4",
  "OUTCOME_5"
)

age_var <- "AGE"
sex_var <- "SEX"
baseline_severity_var <- "BASELINE_SEVERITY"

pc_vars <- paste0("PC", 1:5)

# Site covariates can differ across outcomes.
# Each list entry should contain the site variables retained
# for the corresponding outcome.
site_terms_by_outcome <- list(
  OUTCOME_1 = c("SITE_1", "SITE_2"),
  OUTCOME_2 = c("SITE_1", "SITE_2"),
  OUTCOME_3 = c("SITE_1", "SITE_2"),
  OUTCOME_4 = c("SITE_1", "SITE_2"),
  OUTCOME_5 = c("SITE_1", "SITE_2")
)

delta_aic_cut <- 2
lrt_p_cut <- 0.05

# Out-of-fold AUC settings.
n_folds <- 10
cv_seed <- 123

# Predicted probabilities for retained interactions are
# calculated at -2 and +2 SD on the original PGS scale.
# Because the PGS is Gelman-scaled as (x - mean) / (2 * SD),
# these values correspond to -1 and +1 on the model scale.
prediction_prs_values <- c(-1, 1)

# Continuous moderators are evaluated at the 25th, 50th,
# and 75th percentiles of their raw distributions.
moderator_probs <- c(0.25, 0.50, 0.75)


# ============================================================
# 2. Helper functions
# ============================================================

# Gelman scaling centres continuous predictors at zero and divides
# them by two standard deviations. This makes the scale of continuous
# predictors more comparable with binary predictors while preserving
# their original ordering.
gelman_scale <- function(x) {
  if (all(is.na(x))) return(x)

  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(x - m)
  }

  (x - m) / (2 * s)
}


# Convert raw moderator values (e.g., P25, P50, P75) to the same
# Gelman-scaled metric used when fitting the regression model.
# This is needed when simple slopes are evaluated at meaningful
# values of a continuous moderator.
raw_to_gelman_scaled <- function(x, raw_reference) {
  s <- sd(raw_reference, na.rm = TRUE)
  m <- mean(raw_reference, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(x - m)
  }

  (x - m) / (2 * s)
}


# Check that each outcome is coded as 0/1 before model fitting.
# Character or factor values are converted to numeric when possible;
# any other coding stops the analysis rather than being silently recoded.
as_binary_01 <- function(x, name = "outcome") {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    x <- suppressWarnings(as.numeric(x))
  }

  values <- sort(unique(x[!is.na(x)]))

  if (!all(values %in% c(0, 1))) {
    stop(
      name,
      " must be coded 0/1. Observed values: ",
      paste(values, collapse = ", ")
    )
  }

  as.integer(x)
}


# Build the logistic-regression formula for one PGS-outcome pair.
# The base model always contains the PGS and the prespecified
# covariates. Interaction terms are added only when requested.
build_formula_from_interactions <- function(
  outcome,
  prs,
  interactions,
  site_terms,
  pc_vars,
  age_var,
  sex_var,
  baseline_severity_var
) {

  main_terms <- c(
    prs,
    baseline_severity_var,
    sex_var,
    age_var,
    pc_vars
  )

  rhs_terms <- main_terms

  if ("severity" %in% interactions) {
    rhs_terms <- c(
      rhs_terms,
      paste0(prs, ":", baseline_severity_var)
    )
  }

  if ("sex" %in% interactions) {
    rhs_terms <- c(
      rhs_terms,
      paste0(prs, ":", sex_var)
    )
  }

  if ("age" %in% interactions) {
    rhs_terms <- c(
      rhs_terms,
      paste0(prs, ":", age_var)
    )
  }

  rhs_string <- paste(
    c(rhs_terms, site_terms),
    collapse = " + "
  )

  as.formula(
    paste(outcome, "~", rhs_string)
  )
}


# Use the same complete-case sample for the base model and all
# competing interaction models within a given PGS-outcome pair.
# This avoids comparing AIC or LRT statistics across models fitted
# to different sets of participants.
prepare_complete_case_data <- function(
  dat,
  outcome,
  prs,
  site_terms,
  pc_vars,
  age_var,
  sex_var,
  baseline_severity_var
) {

  needed_vars <- c(
    outcome,
    prs,
    age_var,
    sex_var,
    baseline_severity_var,
    pc_vars,
    site_terms
  )

  dat %>%
    filter(!is.na(.data[[outcome]])) %>%
    mutate(
      "{sex_var}" := factor(.data[[sex_var]])
    ) %>%
    select(all_of(needed_vars)) %>%
    drop_na()
}


extract_lrt_p <- function(anova_object) {
  out <- anova_object$`Pr(>Chi)`[2]

  if (length(out) == 0) {
    return(NA_real_)
  }

  out
}


interaction_label <- function(interactions) {
  if (length(interactions) == 0) {
    return("m_base")
  }

  paste0(
    "m_",
    paste(sort(interactions), collapse = "_")
  )
}


# Retain an interaction only when it satisfies both prespecified
# model-comparison criteria:
#   1. LRT p < lrt_p_cut
#   2. AIC improves by at least delta_aic_cut points
#
# More than one interaction can be retained in the final model if
# each one improves fit relative to the same base model.
select_interactions <- function(
  aic_base,
  aic_severity,
  aic_sex,
  aic_age,
  p_severity,
  p_sex,
  p_age,
  delta_aic_cut = 2,
  lrt_p_cut = 0.05
) {

  keep_severity <- (
    !is.na(p_severity) &&
      p_severity < lrt_p_cut &&
      (aic_base - aic_severity) >= delta_aic_cut
  )

  keep_sex <- (
    !is.na(p_sex) &&
      p_sex < lrt_p_cut &&
      (aic_base - aic_sex) >= delta_aic_cut
  )

  keep_age <- (
    !is.na(p_age) &&
      p_age < lrt_p_cut &&
      (aic_base - aic_age) >= delta_aic_cut
  )

  selected <- c(
    "severity",
    "sex",
    "age"
  )[c(
    keep_severity,
    keep_sex,
    keep_age
  )]

  list(
    selected = selected,
    keep_severity = keep_severity,
    keep_sex = keep_sex,
    keep_age = keep_age
  )
}


# Compute AUC from ranks without requiring an additional ROC package.
# AUC is calculated only when both outcome classes are represented.
compute_auc_manual <- function(y_true, y_score) {
  keep <- !is.na(y_true) & !is.na(y_score)

  y_true <- y_true[keep]
  y_score <- y_score[keep]

  if (length(unique(y_true)) < 2) {
    return(NA_real_)
  }

  y_true <- as.numeric(y_true)

  n1 <- sum(y_true == 1)
  n0 <- sum(y_true == 0)

  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }

  ranks <- rank(
    y_score,
    ties.method = "average"
  )

  (
    sum(ranks[y_true == 1]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}


# Create approximately class-balanced folds so that cases and
# non-cases are distributed across cross-validation folds.
make_stratified_folds <- function(
  y,
  k = 10,
  seed = 123
) {

  set.seed(seed)

  y <- as.numeric(y)

  idx1 <- which(y == 1)
  idx0 <- which(y == 0)

  idx1 <- sample(idx1)
  idx0 <- sample(idx0)

  folds <- integer(length(y))

  folds[idx1] <- rep(
    seq_len(k),
    length.out = length(idx1)
  )

  folds[idx0] <- rep(
    seq_len(k),
    length.out = length(idx0)
  )

  folds
}


# Estimate scaling parameters using the training fold only.
# These values are then applied unchanged to the held-out fold.
fit_scaler <- function(
  train_df,
  scale_vars
) {

  means <- sapply(
    train_df[, scale_vars, drop = FALSE],
    function(x) mean(x, na.rm = TRUE)
  )

  sds <- sapply(
    train_df[, scale_vars, drop = FALSE],
    function(x) sd(x, na.rm = TRUE)
  )

  sds[is.na(sds) | sds == 0] <- NA_real_

  list(
    means = means,
    sds = sds
  )
}


# Apply training-fold scaling parameters to a dataset.
# Keeping scaling fold-specific prevents information from the
# held-out observations from entering model estimation.
apply_scaler <- function(
  df,
  scaler,
  scale_vars
) {

  out <- df

  for (v in scale_vars) {
    mu <- scaler$means[[v]]
    s <- scaler$sds[[v]]

    if (is.na(s)) {
      out[[v]] <- df[[v]] - mu
    } else {
      out[[v]] <- (
        df[[v]] - mu
      ) / (2 * s)
    }
  }

  out
}


# Estimate out-of-fold AUC using stratified k-fold cross-validation.
# Predictions for each participant are generated by a model that was
# fitted without that participant. Scaling is also estimated within
# each training fold to avoid data leakage.
compute_oof_auc <- function(
  formula_object,
  data_raw,
  outcome_name,
  scale_vars,
  sex_var,
  k = 10,
  seed = 123
) {

  y <- data_raw[[outcome_name]]

  if (!all(sort(unique(y)) %in% c(0, 1))) {
    stop(
      "Outcome ",
      outcome_name,
      " must be coded 0/1 for OOF AUC computation."
    )
  }

  folds <- make_stratified_folds(
    y = y,
    k = k,
    seed = seed
  )

  oof_pred <- rep(
    NA_real_,
    nrow(data_raw)
  )

  for (fold in seq_len(k)) {

    test_idx <- which(folds == fold)
    train_idx <- which(folds != fold)

    train_raw <- data_raw[
      train_idx,
      ,
      drop = FALSE
    ]

    test_raw <- data_raw[
      test_idx,
      ,
      drop = FALSE
    ]

    scaler <- fit_scaler(
      train_df = train_raw,
      scale_vars = scale_vars
    )

    train_df <- apply_scaler(
      df = train_raw,
      scaler = scaler,
      scale_vars = scale_vars
    )

    test_df <- apply_scaler(
      df = test_raw,
      scaler = scaler,
      scale_vars = scale_vars
    )

    train_df[[sex_var]] <- factor(
      train_df[[sex_var]]
    )

    test_df[[sex_var]] <- factor(
      test_df[[sex_var]],
      levels = levels(train_df[[sex_var]])
    )

    fit <- glm(
      formula_object,
      family = binomial(),
      data = train_df
    )

    oof_pred[test_idx] <- predict(
      fit,
      newdata = test_df,
      type = "response"
    )
  }

  auc <- compute_auc_manual(
    y_true = y,
    y_score = oof_pred
  )

  list(
    auc = auc,
    oof_pred = oof_pred,
    folds = folds
  )
}


# Define a common reference profile for model-based probabilities.
# Continuous covariates are fixed at their medians, sex at the
# reference factor level, and dummy-coded site variables at zero.
get_reference_profile <- function(
  dat,
  age_var,
  sex_var,
  baseline_severity_var,
  pc_vars,
  site_terms
) {

  ref <- list()

  ref[[age_var]] <- median(
    dat[[age_var]],
    na.rm = TRUE
  )

  ref[[baseline_severity_var]] <- median(
    dat[[baseline_severity_var]],
    na.rm = TRUE
  )

  ref[[sex_var]] <- levels(
    dat[[sex_var]]
  )[1]

  for (pc in pc_vars) {
    ref[[pc]] <- median(
      dat[[pc]],
      na.rm = TRUE
    )
  }

  for (site in site_terms) {
    ref[[site]] <- 0
  }

  ref
}


# Construct post-hoc prediction grids for retained interactions.
# For continuous moderators, predictions are evaluated at P25,
# P50, and P75. For sex, predictions are evaluated within each
# observed factor level. The PGS is fixed at -2 and +2 raw SD,
# corresponding to -1 and +1 after Gelman scaling.
build_prediction_data <- function(
  dat,
  prs,
  interaction,
  site_terms,
  pc_vars,
  age_var,
  sex_var,
  baseline_severity_var,
  prs_values = c(-1, 1)
) {

  ref <- get_reference_profile(
    dat = dat,
    age_var = age_var,
    sex_var = sex_var,
    baseline_severity_var = baseline_severity_var,
    pc_vars = pc_vars,
    site_terms = site_terms
  )

  if (interaction == "sex") {

    nd <- expand.grid(
      PRS_VALUE = prs_values,
      SEX_LEVEL = levels(dat[[sex_var]]),
      stringsAsFactors = FALSE
    )

    nd[[prs]] <- nd$PRS_VALUE

    nd[[age_var]] <- ref[[age_var]]
    nd[[baseline_severity_var]] <- ref[[baseline_severity_var]]

    nd[[sex_var]] <- factor(
      nd$SEX_LEVEL,
      levels = levels(dat[[sex_var]])
    )

    nd$moderator_level <- paste0(
      sex_var,
      " = ",
      nd$SEX_LEVEL
    )
  }

  if (interaction == "age") {

    moderator_values <- as.numeric(
      quantile(
        dat[[age_var]],
        probs = moderator_probs,
        na.rm = TRUE
      )
    )

    moderator_labels <- c(
      "P25",
      "P50",
      "P75"
    )

    nd <- expand.grid(
      PRS_VALUE = prs_values,
      moderator_level = moderator_labels,
      stringsAsFactors = FALSE
    )

    nd[[prs]] <- nd$PRS_VALUE

    nd[[age_var]] <- moderator_values[
      match(
        nd$moderator_level,
        moderator_labels
      )
    ]

    nd[[baseline_severity_var]] <- ref[[baseline_severity_var]]

    nd[[sex_var]] <- factor(
      ref[[sex_var]],
      levels = levels(dat[[sex_var]])
    )
  }

  if (interaction == "severity") {

    moderator_values <- as.numeric(
      quantile(
        dat[[baseline_severity_var]],
        probs = moderator_probs,
        na.rm = TRUE
      )
    )

    moderator_labels <- c(
      "P25",
      "P50",
      "P75"
    )

    nd <- expand.grid(
      PRS_VALUE = prs_values,
      moderator_level = moderator_labels,
      stringsAsFactors = FALSE
    )

    nd[[prs]] <- nd$PRS_VALUE

    nd[[age_var]] <- ref[[age_var]]

    nd[[baseline_severity_var]] <- moderator_values[
      match(
        nd$moderator_level,
        moderator_labels
      )
    ]

    nd[[sex_var]] <- factor(
      ref[[sex_var]],
      levels = levels(dat[[sex_var]])
    )
  }

  for (pc in pc_vars) {
    nd[[pc]] <- ref[[pc]]
  }

  for (site in site_terms) {
    nd[[site]] <- ref[[site]]
  }

  nd
}


# Obtain model-based predicted probabilities and Wald 95% confidence
# intervals on the probability scale. These values are used to
# describe fitted interactions rather than as individual predictions.
add_prediction_intervals <- function(
  fit,
  newdata
) {

  pred <- predict(
    fit,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )

  newdata %>%
    mutate(
      eta = pred$fit,
      eta_se = pred$se.fit,
      predicted_prob = plogis(eta),
      CI_low_prob = plogis(
        eta - 1.96 * eta_se
      ),
      CI_high_prob = plogis(
        eta + 1.96 * eta_se
      )
    )
}


# Format emtrends output for simple-slope analyses.
# The PGS slope is estimated on the log-odds scale and exponentiated
# to obtain an odds ratio at each moderator level.
clean_emtrends_summary <- function(
  x,
  prs,
  outcome,
  moderator,
  model_n,
  raw_levels = NULL
) {

  out <- as.data.frame(x)

  trend_col <- paste0(
    prs,
    ".trend"
  )

  names(out)[names(out) == trend_col] <- "logit_slope"
  names(out)[names(out) == "asymp.LCL"] <- "logit_CI_low"
  names(out)[names(out) == "asymp.UCL"] <- "logit_CI_high"
  names(out)[names(out) == "z.ratio"] <- "z_ratio"
  names(out)[names(out) == "p.value"] <- "p_value"

  out <- out %>%
    mutate(
      outcome = outcome,
      PGS = prs,
      moderator = moderator,
      n_complete_cases = model_n,
      OR = exp(logit_slope),
      OR_CI_low = exp(logit_CI_low),
      OR_CI_high = exp(logit_CI_high)
    )

  if (!is.null(raw_levels)) {

    out <- out %>%
      mutate(
        moderator_level = raw_levels$level_label[
          seq_len(n())
        ],
        moderator_raw_value = raw_levels$raw_value[
          seq_len(n())
        ],
        moderator_scaled_value = raw_levels$scaled_value[
          seq_len(n())
        ]
      )

  } else {

    moderator_col <- intersect(
      names(out),
      moderator
    )

    if (length(moderator_col) == 1) {
      out$moderator_level <- paste0(
        moderator,
        " = ",
        out[[moderator_col]]
      )
    } else {
      out$moderator_level <- NA_character_
    }

    out$moderator_raw_value <- NA_real_
    out$moderator_scaled_value <- NA_real_
  }

  out %>%
    select(
      outcome,
      PGS,
      moderator,
      moderator_level,
      moderator_raw_value,
      moderator_scaled_value,
      n_complete_cases,
      logit_slope,
      SE,
      df,
      logit_CI_low,
      logit_CI_high,
      z_ratio,
      p_value,
      OR,
      OR_CI_low,
      OR_CI_high,
      everything()
    )
}


# Format pairwise contrasts between simple PGS slopes.
# Exponentiating the contrast gives the ratio between the PGS odds
# ratios estimated at two moderator levels.
clean_contrast_summary <- function(
  x,
  prs,
  outcome,
  moderator,
  model_n,
  contrast_labels = NULL
) {

  out <- as.data.frame(x)

  names(out)[names(out) == "asymp.LCL"] <- "log_OR_ratio_CI_low"
  names(out)[names(out) == "asymp.UCL"] <- "log_OR_ratio_CI_high"
  names(out)[names(out) == "z.ratio"] <- "z_ratio"
  names(out)[names(out) == "p.value"] <- "p_value"
  names(out)[names(out) == "estimate"] <- "log_OR_ratio"

  if (
    !is.null(contrast_labels) &&
      length(contrast_labels) == nrow(out)
  ) {
    out$contrast_readable <- contrast_labels
  } else {
    out$contrast_readable <- out$contrast
  }

  out %>%
    mutate(
      outcome = outcome,
      PGS = prs,
      moderator = moderator,
      n_complete_cases = model_n,
      OR_ratio = exp(log_OR_ratio),
      OR_ratio_CI_low = exp(log_OR_ratio_CI_low),
      OR_ratio_CI_high = exp(log_OR_ratio_CI_high)
    ) %>%
    select(
      outcome,
      PGS,
      moderator,
      contrast_readable,
      contrast,
      n_complete_cases,
      log_OR_ratio,
      SE,
      df,
      log_OR_ratio_CI_low,
      log_OR_ratio_CI_high,
      z_ratio,
      p_value,
      OR_ratio,
      OR_ratio_CI_low,
      OR_ratio_CI_high,
      everything()
    )
}


# ============================================================
# 3. Read and prepare data
# ============================================================

# Keep one unscaled copy for fold-specific cross-validation and for
# recovering raw moderator percentiles, and one globally scaled copy
# for coefficient estimation and post-hoc analyses in the final models.
dat_raw <- read.csv(
  data_file,
  stringsAsFactors = FALSE
)

# Optional generic validity filter.
# Modify or remove if age values are handled differently.
dat_raw <- dat_raw %>%
  filter(
    is.na(.data[[age_var]]) |
      .data[[age_var]] > 0
  )

for (outcome in outcomes) {
  dat_raw[[outcome]] <- as_binary_01(
    dat_raw[[outcome]],
    name = outcome
  )
}

continuous_vars <- c(
  age_var,
  baseline_severity_var,
  prs_vars,
  pc_vars
)

existing_continuous_vars <- intersect(
  continuous_vars,
  names(dat_raw)
)

missing_prs_vars <- setdiff(
  prs_vars,
  names(dat_raw)
)

if (length(missing_prs_vars) > 0) {
  stop(
    "Missing PGS variables: ",
    paste(
      missing_prs_vars,
      collapse = ", "
    )
  )
}

dat_scaled <- dat_raw %>%
  mutate(
    across(
      all_of(existing_continuous_vars),
      gelman_scale
    )
  )


# ============================================================
# 4. Fit candidate models and select final models
# ============================================================

# For each outcome and PGS, four candidate models are fitted on the
# same complete-case sample:
#   - base model: PGS + covariates
#   - PGS x baseline severity
#   - PGS x sex
#   - PGS x age
#
# Each interaction model is compared with the base model using LRT
# and AIC. Interactions meeting both criteria are carried forward
# into a single final model.
model_selection_rows <- list()
final_model_rows <- list()
final_models <- list()
oof_prediction_rows <- list()

for (outcome in outcomes) {

  site_terms <- site_terms_by_outcome[[outcome]]

  if (is.null(site_terms)) {
    stop(
      "No site-variable specification found for ",
      outcome
    )
  }

  for (prs in prs_vars) {

    dat_cc <- prepare_complete_case_data(
      dat = dat_scaled,
      outcome = outcome,
      prs = prs,
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    dat_cc_raw <- prepare_complete_case_data(
      dat = dat_raw,
      outcome = outcome,
      prs = prs,
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    n_cc <- nrow(dat_cc)

    if (n_cc == 0) {
      warning(
        "No complete-case data for ",
        outcome,
        " and ",
        prs
      )
      next
    }

    f_base <- build_formula_from_interactions(
      outcome = outcome,
      prs = prs,
      interactions = character(0),
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    f_severity <- build_formula_from_interactions(
      outcome = outcome,
      prs = prs,
      interactions = "severity",
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    f_sex <- build_formula_from_interactions(
      outcome = outcome,
      prs = prs,
      interactions = "sex",
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    f_age <- build_formula_from_interactions(
      outcome = outcome,
      prs = prs,
      interactions = "age",
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    m_base <- glm(
      f_base,
      family = binomial(),
      data = dat_cc
    )

    m_severity <- glm(
      f_severity,
      family = binomial(),
      data = dat_cc
    )

    m_sex <- glm(
      f_sex,
      family = binomial(),
      data = dat_cc
    )

    m_age <- glm(
      f_age,
      family = binomial(),
      data = dat_cc
    )

    aic_base <- AIC(m_base)
    aic_severity <- AIC(m_severity)
    aic_sex <- AIC(m_sex)
    aic_age <- AIC(m_age)

    lrt_severity <- anova(
      m_base,
      m_severity,
      test = "LRT"
    )

    lrt_sex <- anova(
      m_base,
      m_sex,
      test = "LRT"
    )

    lrt_age <- anova(
      m_base,
      m_age,
      test = "LRT"
    )

    p_severity <- extract_lrt_p(
      lrt_severity
    )

    p_sex <- extract_lrt_p(
      lrt_sex
    )

    p_age <- extract_lrt_p(
      lrt_age
    )

    selected <- select_interactions(
      aic_base = aic_base,
      aic_severity = aic_severity,
      aic_sex = aic_sex,
      aic_age = aic_age,
      p_severity = p_severity,
      p_sex = p_sex,
      p_age = p_age,
      delta_aic_cut = delta_aic_cut,
      lrt_p_cut = lrt_p_cut
    )

    final_interactions <- selected$selected

    chosen_model <- interaction_label(
      final_interactions
    )

    f_final <- build_formula_from_interactions(
      outcome = outcome,
      prs = prs,
      interactions = final_interactions,
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    m_final <- glm(
      f_final,
      family = binomial(),
      data = dat_cc
    )

    # --------------------------------------------------------
    # Out-of-fold AUC
    # --------------------------------------------------------
    #
    # OOF AUC is reported for the base and selected final model.
    # This measures discrimination of the full multivariable models,
    # not of the PGS term in isolation. The AUC difference is kept as
    # a descriptive measure of incremental discrimination.

    scale_vars_cv <- intersect(
      c(
        age_var,
        baseline_severity_var,
        prs,
        pc_vars
      ),
      names(dat_cc_raw)
    )

    cv_base <- compute_oof_auc(
      formula_object = f_base,
      data_raw = dat_cc_raw,
      outcome_name = outcome,
      scale_vars = scale_vars_cv,
      sex_var = sex_var,
      k = n_folds,
      seed = cv_seed
    )

    cv_final <- compute_oof_auc(
      formula_object = f_final,
      data_raw = dat_cc_raw,
      outcome_name = outcome,
      scale_vars = scale_vars_cv,
      sex_var = sex_var,
      k = n_folds,
      seed = cv_seed
    )

    auc_oof_base <- cv_base$auc
    auc_oof_final <- cv_final$auc
    delta_auc_oof <- (
      auc_oof_final -
        auc_oof_base
    )

    # --------------------------------------------------------
    # Final model coefficient table
    # --------------------------------------------------------
    #
    # Wald confidence intervals are used so that coefficient, OR,
    # CI, and p-value estimates are obtained consistently for all
    # terms in the selected frequentist model.

    coef_tab <- as.data.frame(
      summary(m_final)$coefficients
    )

    coef_tab$term <- rownames(coef_tab)
    rownames(coef_tab) <- NULL

    names(coef_tab) <- c(
      "estimate",
      "SE",
      "z_value",
      "p_value",
      "term"
    )

    ci_mat <- confint.default(
      m_final
    )

    ci_tab <- data.frame(
      term = rownames(ci_mat),
      CI_low_beta = ci_mat[, 1],
      CI_high_beta = ci_mat[, 2],
      stringsAsFactors = FALSE
    )

    final_tab <- coef_tab %>%
      left_join(
        ci_tab,
        by = "term"
      ) %>%
      mutate(
        OR = exp(estimate),
        CI_low = exp(CI_low_beta),
        CI_high = exp(CI_high_beta),
        outcome = outcome,
        PGS = prs,
        chosen_model = chosen_model,
        n_complete_cases = n_cc,
        AUC_OOF_base = auc_oof_base,
        AUC_OOF_final = auc_oof_final,
        delta_AUC_OOF = delta_auc_oof
      ) %>%
      select(
        outcome,
        PGS,
        chosen_model,
        n_complete_cases,
        AUC_OOF_base,
        AUC_OOF_final,
        delta_AUC_OOF,
        term,
        estimate,
        SE,
        z_value,
        p_value,
        OR,
        CI_low,
        CI_high
      )

    model_selection_rows[[
      paste(
        outcome,
        prs,
        sep = "__"
      )
    ]] <- tibble(
      outcome = outcome,
      PGS = prs,
      n_complete_cases = n_cc,
      site_terms = paste(
        site_terms,
        collapse = " + "
      ),
      AIC_base = aic_base,
      AIC_severity = aic_severity,
      AIC_sex = aic_sex,
      AIC_age = aic_age,
      p_LRT_severity = p_severity,
      p_LRT_sex = p_sex,
      p_LRT_age = p_age,
      deltaAIC_severity = (
        aic_base -
          aic_severity
      ),
      deltaAIC_sex = (
        aic_base -
          aic_sex
      ),
      deltaAIC_age = (
        aic_base -
          aic_age
      ),
      keep_severity = selected$keep_severity,
      keep_sex = selected$keep_sex,
      keep_age = selected$keep_age,
      chosen_model = chosen_model,
      selected_interactions = ifelse(
        length(final_interactions) == 0,
        "none",
        paste(
          final_interactions,
          collapse = ","
        )
      ),
      AUC_OOF_base = auc_oof_base,
      AUC_OOF_final = auc_oof_final,
      delta_AUC_OOF = delta_auc_oof
    )

    final_model_rows[[
      paste(
        outcome,
        prs,
        sep = "__"
      )
    ]] <- final_tab

    final_models[[
      paste(
        outcome,
        prs,
        sep = "__"
      )
    ]] <- list(
      fit = m_final,
      data = dat_cc,
      data_raw = dat_cc_raw,
      outcome = outcome,
      prs = prs,
      chosen_model = chosen_model,
      selected_interactions = final_interactions,
      site_terms = site_terms
    )

    oof_prediction_rows[[
      paste(
        outcome,
        prs,
        "base",
        sep = "__"
      )
    ]] <- tibble(
      outcome = outcome,
      PGS = prs,
      model = "base",
      y_true = dat_cc_raw[[outcome]],
      oof_pred = cv_base$oof_pred
    )

    oof_prediction_rows[[
      paste(
        outcome,
        prs,
        "final",
        sep = "__"
      )
    ]] <- tibble(
      outcome = outcome,
      PGS = prs,
      model = "final",
      y_true = dat_cc_raw[[outcome]],
      oof_pred = cv_final$oof_pred
    )
  }
}

model_selection_summary <- bind_rows(
  model_selection_rows
)

final_OR_table <- bind_rows(
  final_model_rows
)

oof_predictions_table <- bind_rows(
  oof_prediction_rows
)


# ============================================================
# 5. Retained interaction models
# ============================================================

# Only models containing interactions selected by the LRT + AIC
# procedure are taken forward to the post-hoc interaction analyses.
selected_interaction_models <- model_selection_summary %>%
  filter(
    selected_interactions != "none"
  )


# ============================================================
# 6. Predicted probabilities for retained interactions
# ============================================================

# Predicted probabilities are used to display the fitted interaction
# pattern on an interpretable probability scale.
#
# Continuous moderators:
#   P25, P50, and P75 of the moderator
#
# PGS:
#   -2 and +2 SD on the original PGS scale
#
# Other covariates are held at the reference profile defined above.
predicted_probability_rows <- list()

for (
  i in seq_len(
    nrow(selected_interaction_models)
  )
) {

  outcome <- selected_interaction_models$outcome[i]
  prs <- selected_interaction_models$PGS[i]

  interactions <- unlist(
    strsplit(
      selected_interaction_models$selected_interactions[i],
      "\\s*,\\s*"
    )
  )

  key <- paste(
    outcome,
    prs,
    sep = "__"
  )

  obj <- final_models[[key]]

  for (interaction in interactions) {

    nd <- build_prediction_data(
      dat = obj$data,
      prs = prs,
      interaction = interaction,
      site_terms = obj$site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var,
      prs_values = prediction_prs_values
    )

    nd <- add_prediction_intervals(
      fit = obj$fit,
      newdata = nd
    )

    nd <- nd %>%
      mutate(
        outcome = outcome,
        PGS = prs,
        chosen_model = obj$chosen_model,
        interaction = interaction
      ) %>%
      select(
        outcome,
        PGS,
        chosen_model,
        interaction,
        moderator_level,
        PRS_VALUE,
        predicted_prob,
        CI_low_prob,
        CI_high_prob,
        everything()
      )

    predicted_probability_rows[[
      paste(
        outcome,
        prs,
        interaction,
        sep = "__"
      )
    ]] <- nd
  }
}

predicted_probabilities <- bind_rows(
  predicted_probability_rows
)


# ============================================================
# 7. Simple-slope analyses and contrasts
# ============================================================

# emtrends() is used to estimate the conditional PGS slope from each
# retained interaction model.
#
# For continuous moderators, the PGS slope is estimated at P25,
# P50, and P75. For sex interactions, the slope is estimated within
# each sex category.
#
# The simple-slope p-value tests whether the PGS-outcome association
# differs from zero at that moderator level.
#
# Pairwise contrasts then test whether the PGS slopes differ across
# moderator levels. No additional p-value adjustment is applied to
# these descriptive post-hoc contrasts.
simple_slope_rows <- list()
slope_contrast_rows <- list()
interaction_model_audit_rows <- list()

for (
  i in seq_len(
    nrow(selected_interaction_models)
  )
) {

  outcome <- selected_interaction_models$outcome[i]
  prs <- selected_interaction_models$PGS[i]

  interactions <- unlist(
    strsplit(
      selected_interaction_models$selected_interactions[i],
      "\\s*,\\s*"
    )
  )

  key <- paste(
    outcome,
    prs,
    sep = "__"
  )

  obj <- final_models[[key]]

  fit <- obj$fit
  dat_cc <- obj$data
  dat_cc_raw <- obj$data_raw

  interaction_model_audit_rows[[
    length(
      interaction_model_audit_rows
    ) + 1
  ]] <- tibble(
    outcome = outcome,
    PGS = prs,
    selected_interactions = paste(
      interactions,
      collapse = ","
    ),
    n_complete_cases = nrow(dat_cc),
    formula = paste(
      deparse(
        formula(fit)
      ),
      collapse = " "
    )
  )

  for (interaction in interactions) {

    # --------------------------------------------------------
    # Sex interaction
    # --------------------------------------------------------

    if (interaction == "sex") {

      specs_formula <- as.formula(
        paste(
          "~",
          sex_var
        )
      )

      trends <- emtrends(
        fit,
        specs = specs_formula,
        var = prs
      )

      trend_summary <- summary(
        trends,
        infer = c(TRUE, TRUE),
        type = "link"
      )

      simple_slope_rows[[
        length(simple_slope_rows) + 1
      ]] <- clean_emtrends_summary(
        x = trend_summary,
        prs = prs,
        outcome = outcome,
        moderator = sex_var,
        model_n = nrow(dat_cc)
      )

      contrasts <- pairs(
        trends,
        adjust = "none"
      )

      contrast_summary <- summary(
        contrasts,
        infer = c(TRUE, TRUE),
        type = "link"
      )

      slope_contrast_rows[[
        length(slope_contrast_rows) + 1
      ]] <- clean_contrast_summary(
        x = contrast_summary,
        prs = prs,
        outcome = outcome,
        moderator = sex_var,
        model_n = nrow(dat_cc)
      )
    }

    # --------------------------------------------------------
    # Continuous-moderator interactions
    # --------------------------------------------------------

    if (
      interaction %in%
        c(
          "age",
          "severity"
        )
    ) {

      moderator <- ifelse(
        interaction == "age",
        age_var,
        baseline_severity_var
      )

      raw_values <- as.numeric(
        quantile(
          dat_cc_raw[[moderator]],
          probs = moderator_probs,
          na.rm = TRUE
        )
      )

      scaled_values <- raw_to_gelman_scaled(
        x = raw_values,
        raw_reference = dat_cc_raw[[moderator]]
      )

      raw_levels <- tibble(
        level_label = c(
          "P25",
          "P50",
          "P75"
        ),
        raw_value = raw_values,
        scaled_value = scaled_values
      )

      at_list <- list(
        scaled_values
      )

      names(at_list) <- moderator

      specs_formula <- as.formula(
        paste(
          "~",
          moderator
        )
      )

      trends <- emtrends(
        fit,
        specs = specs_formula,
        var = prs,
        at = at_list
      )

      trend_summary <- summary(
        trends,
        infer = c(TRUE, TRUE),
        type = "link"
      )

      simple_slope_rows[[
        length(simple_slope_rows) + 1
      ]] <- clean_emtrends_summary(
        x = trend_summary,
        prs = prs,
        outcome = outcome,
        moderator = moderator,
        model_n = nrow(dat_cc),
        raw_levels = raw_levels
      )

      contrasts <- pairs(
        trends,
        adjust = "none"
      )

      contrast_summary <- summary(
        contrasts,
        infer = c(TRUE, TRUE),
        type = "link"
      )

      contrast_labels <- c(
        paste0(
          "P25 vs P50; ",
          moderator,
          " = ",
          raw_values[1],
          " vs ",
          raw_values[2]
        ),
        paste0(
          "P25 vs P75; ",
          moderator,
          " = ",
          raw_values[1],
          " vs ",
          raw_values[3]
        ),
        paste0(
          "P50 vs P75; ",
          moderator,
          " = ",
          raw_values[2],
          " vs ",
          raw_values[3]
        )
      )

      slope_contrast_rows[[
        length(slope_contrast_rows) + 1
      ]] <- clean_contrast_summary(
        x = contrast_summary,
        prs = prs,
        outcome = outcome,
        moderator = moderator,
        model_n = nrow(dat_cc),
        contrast_labels = contrast_labels
      )
    }
  }
}

simple_prs_slopes <- bind_rows(
  simple_slope_rows
)

prs_slope_contrasts <- bind_rows(
  slope_contrast_rows
)

interaction_model_audit <- bind_rows(
  interaction_model_audit_rows
)


# ============================================================
# 8. Save outputs
# ============================================================

# Results are written both as separate CSV files and as a single
# Excel workbook with one worksheet per analysis component.
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  model_selection_summary,
  file.path(
    output_dir,
    "model_selection_summary.csv"
  )
)

write_csv(
  final_OR_table,
  file.path(
    output_dir,
    "final_OR_table.csv"
  )
)

write_csv(
  oof_predictions_table,
  file.path(
    output_dir,
    "oof_predictions_table.csv"
  )
)

write_csv(
  selected_interaction_models,
  file.path(
    output_dir,
    "selected_interaction_models.csv"
  )
)

write_csv(
  predicted_probabilities,
  file.path(
    output_dir,
    "predicted_probabilities.csv"
  )
)

write_csv(
  simple_prs_slopes,
  file.path(
    output_dir,
    "simple_prs_slopes.csv"
  )
)

write_csv(
  prs_slope_contrasts,
  file.path(
    output_dir,
    "prs_slope_contrasts.csv"
  )
)

write_csv(
  interaction_model_audit,
  file.path(
    output_dir,
    "interaction_model_audit.csv"
  )
)

write_xlsx(
  list(
    model_selection = model_selection_summary,
    final_models = final_OR_table,
    OOF_predictions = oof_predictions_table,
    selected_interactions = selected_interaction_models,
    predicted_probabilities = predicted_probabilities,
    simple_slopes = simple_prs_slopes,
    slope_contrasts = prs_slope_contrasts,
    interaction_model_audit = interaction_model_audit
  ),
  path = file.path(
    output_dir,
    "frequentist_analysis_outputs.xlsx"
  )
)

cat(
  "\nFrequentist analysis outputs saved in:\n",
  output_dir,
  "\n"
)
