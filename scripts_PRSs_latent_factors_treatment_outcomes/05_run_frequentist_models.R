# ============================================================
# Project: Latent Genomic Factors and Treatment Outcomes
# Script: 05_run_frequentist_models.R
#
# Purpose:
# - fit logistic regression models for multiple binary outcomes
# - test associations between polygenic scores (PGSs) and outcomes
# - evaluate PGS-by-moderator interaction terms
# - compare nested models using likelihood-ratio tests and AIC
# - estimate out-of-fold AUC with fold-specific scaling
# - quantify incremental model performance using Tjur's and
#   McFadden's pseudo-R2 relative to a covariate-only model
# - calculate predicted probabilities for significant retained
#   interactions
# - estimate simple PGS slopes and contrasts between slopes
# - save analysis-ready summary tables
#
# Notes:
# - this script is a generic template and should be adapted to
#   the variable names and structure of the cohort
# - continuous variables are scaled using Gelman scaling:
#   (x - mean) / (2 * SD)
# - each PGS-by-moderator interaction model is evaluated
#   separately against the base model
# - an interaction term is retained when both:
#     LRT p < 0.05
#     delta AIC >= 2 relative to the base model
# - if more than one interaction meets both criteria, all
#   qualifying interaction terms are included jointly in the
#   final model
# - post-hoc analyses are performed for retained interaction
#   terms with Wald p < 0.05 in the final model
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
pgs_vars <- c(
  "PGS_1",
  "PGS_2",
  "PGS_3",
  "PGS_4",
  "PGS_5"
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
prediction_pgs_values <- c(-1, 1)

# Continuous moderators are evaluated at the 25th, 50th,
# and 75th percentiles of their raw distributions.
moderator_probs <- c(0.25, 0.50, 0.75)


# ============================================================
# 2. Helper functions
# ============================================================

gelman_scale <- function(x) {
  if (all(is.na(x))) return(x)

  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(x - m)
  }

  (x - m) / (2 * s)
}


raw_to_gelman_scaled <- function(x, raw_reference) {
  s <- sd(raw_reference, na.rm = TRUE)
  m <- mean(raw_reference, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(x - m)
  }

  (x - m) / (2 * s)
}


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


build_covariate_formula <- function(
  outcome,
  site_terms,
  pc_vars,
  age_var,
  sex_var,
  baseline_severity_var
) {

  # This model contains the prespecified non-PGS covariates only.
  # It serves as the reference model for incremental pseudo-R2.
  rhs_terms <- c(
    baseline_severity_var,
    sex_var,
    age_var,
    pc_vars,
    site_terms
  )

  as.formula(
    paste(
      outcome,
      "~",
      paste(rhs_terms, collapse = " + ")
    )
  )
}


build_formula_from_interactions <- function(
  outcome,
  pgs,
  interactions,
  site_terms,
  pc_vars,
  age_var,
  sex_var,
  baseline_severity_var
) {

  main_terms <- c(
    pgs,
    baseline_severity_var,
    sex_var,
    age_var,
    pc_vars
  )

  rhs_terms <- main_terms

  if ("severity" %in% interactions) {
    rhs_terms <- c(
      rhs_terms,
      paste0(pgs, ":", baseline_severity_var)
    )
  }

  if ("sex" %in% interactions) {
    rhs_terms <- c(
      rhs_terms,
      paste0(pgs, ":", sex_var)
    )
  }

  if ("age" %in% interactions) {
    rhs_terms <- c(
      rhs_terms,
      paste0(pgs, ":", age_var)
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


prepare_complete_case_data <- function(
  dat,
  outcome,
  pgs,
  site_terms,
  pc_vars,
  age_var,
  sex_var,
  baseline_severity_var
) {

  needed_vars <- c(
    outcome,
    pgs,
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


identify_interaction_type <- function(
  term,
  pgs,
  age_var,
  sex_var,
  baseline_severity_var
) {

  age_terms <- c(
    paste0(pgs, ":", age_var),
    paste0(age_var, ":", pgs)
  )

  severity_terms <- c(
    paste0(pgs, ":", baseline_severity_var),
    paste0(baseline_severity_var, ":", pgs)
  )

  sex_prefixes <- c(
    paste0(pgs, ":", sex_var),
    paste0(sex_var, ":", pgs)
  )

  if (term %in% age_terms) {
    return("age")
  }

  if (term %in% severity_terms) {
    return("severity")
  }

  if (any(startsWith(term, sex_prefixes))) {
    return("sex")
  }

  NA_character_
}


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


compute_tjur_r2 <- function(
  fit,
  data,
  outcome_name
) {

  # Tjur's coefficient of discrimination is the difference in
  # mean predicted probability between observed cases and non-cases.
  predicted_prob <- predict(
    fit,
    newdata = data,
    type = "response"
  )

  y <- data[[outcome_name]]

  if (length(unique(y)) < 2) {
    return(NA_real_)
  }

  mean(
    predicted_prob[y == 1],
    na.rm = TRUE
  ) -
    mean(
      predicted_prob[y == 0],
      na.rm = TRUE
    )
}


compute_mcfadden_r2 <- function(fit) {

  # McFadden's pseudo-R2 compares the fitted model log-likelihood
  # with that of an intercept-only model fitted to the same sample.
  null_fit <- update(
    fit,
    formula = . ~ 1
  )

  ll_model <- as.numeric(
    logLik(fit)
  )

  ll_null <- as.numeric(
    logLik(null_fit)
  )

  if (is.na(ll_null) || ll_null == 0) {
    return(NA_real_)
  }

  1 - (ll_model / ll_null)
}


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


build_prediction_data <- function(
  dat,
  pgs,
  interaction,
  site_terms,
  pc_vars,
  age_var,
  sex_var,
  baseline_severity_var,
  pgs_values = c(-1, 1)
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
      PGS_VALUE = pgs_values,
      SEX_LEVEL = levels(dat[[sex_var]]),
      stringsAsFactors = FALSE
    )

    nd[[pgs]] <- nd$PGS_VALUE

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
      PGS_VALUE = pgs_values,
      moderator_level = moderator_labels,
      stringsAsFactors = FALSE
    )

    nd[[pgs]] <- nd$PGS_VALUE

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
      PGS_VALUE = pgs_values,
      moderator_level = moderator_labels,
      stringsAsFactors = FALSE
    )

    nd[[pgs]] <- nd$PGS_VALUE

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


clean_emtrends_summary <- function(
  x,
  pgs,
  outcome,
  moderator,
  model_n,
  raw_levels = NULL
) {

  out <- as.data.frame(x)

  trend_col <- paste0(
    pgs,
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
      PGS = pgs,
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


clean_contrast_summary <- function(
  x,
  pgs,
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
      PGS = pgs,
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
  pgs_vars,
  pc_vars
)

existing_continuous_vars <- intersect(
  continuous_vars,
  names(dat_raw)
)

missing_pgs_vars <- setdiff(
  pgs_vars,
  names(dat_raw)
)

if (length(missing_pgs_vars) > 0) {
  stop(
    "Missing PGS variables: ",
    paste(
      missing_pgs_vars,
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

model_selection_rows <- list()
final_model_rows <- list()
model_performance_rows <- list()
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

  for (pgs in pgs_vars) {

    dat_cc <- prepare_complete_case_data(
      dat = dat_scaled,
      outcome = outcome,
      pgs = pgs,
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    dat_cc_raw <- prepare_complete_case_data(
      dat = dat_raw,
      outcome = outcome,
      pgs = pgs,
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
        pgs
      )
      next
    }

    # --------------------------------------------------------
    # Covariate-only reference model
    # --------------------------------------------------------

    # The covariate-only model is fitted on exactly the same
    # complete-case sample as the candidate and selected models.
    # It provides the reference for the incremental contribution
    # of PGS-related terms to Tjur's and McFadden's pseudo-R2.
    f_covariates <- build_covariate_formula(
      outcome = outcome,
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    f_base <- build_formula_from_interactions(
      outcome = outcome,
      pgs = pgs,
      interactions = character(0),
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    f_severity <- build_formula_from_interactions(
      outcome = outcome,
      pgs = pgs,
      interactions = "severity",
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    f_sex <- build_formula_from_interactions(
      outcome = outcome,
      pgs = pgs,
      interactions = "sex",
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    f_age <- build_formula_from_interactions(
      outcome = outcome,
      pgs = pgs,
      interactions = "age",
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    m_covariates <- glm(
      f_covariates,
      family = binomial(),
      data = dat_cc
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

    # Each interaction model is evaluated separately against
    # the same base model. All interactions meeting both
    # prespecified selection criteria are retained and then
    # included jointly in the final model.
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
      pgs = pgs,
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
    # Incremental pseudo-R2 relative to covariate-only model
    # --------------------------------------------------------

    # For a selected main-effect model, delta pseudo-R2 reflects
    # the addition of the PGS term to the covariate-only model.
    # For a selected moderation model, it reflects the joint
    # addition of the PGS main effect and retained PGS x moderator
    # interaction term(s). Values multiplied by 100 are reported
    # in percentage points (pp), matching the manuscript table.
    tjur_r2_covariates <- compute_tjur_r2(
      fit = m_covariates,
      data = dat_cc,
      outcome_name = outcome
    )

    tjur_r2_selected <- compute_tjur_r2(
      fit = m_final,
      data = dat_cc,
      outcome_name = outcome
    )

    mcfadden_r2_covariates <- compute_mcfadden_r2(
      m_covariates
    )

    mcfadden_r2_selected <- compute_mcfadden_r2(
      m_final
    )

    delta_tjur_r2 <- (
      tjur_r2_selected -
        tjur_r2_covariates
    )

    delta_mcfadden_r2 <- (
      mcfadden_r2_selected -
        mcfadden_r2_covariates
    )

    delta_tjur_r2_pp <- 100 * delta_tjur_r2
    delta_mcfadden_r2_pp <- 100 * delta_mcfadden_r2

    # --------------------------------------------------------
    # Out-of-fold AUC
    # --------------------------------------------------------

    scale_vars_cv <- intersect(
      c(
        age_var,
        baseline_severity_var,
        pgs,
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
        PGS = pgs,
        chosen_model = chosen_model,
        n_complete_cases = n_cc,
        AUC_OOF_base = auc_oof_base,
        AUC_OOF_final = auc_oof_final,
        delta_AUC_OOF = delta_auc_oof,
        Tjur_R2_covariates = tjur_r2_covariates,
        Tjur_R2_selected = tjur_r2_selected,
        delta_Tjur_R2 = delta_tjur_r2,
        delta_Tjur_R2_pp = delta_tjur_r2_pp,
        McFadden_R2_covariates = mcfadden_r2_covariates,
        McFadden_R2_selected = mcfadden_r2_selected,
        delta_McFadden_R2 = delta_mcfadden_r2,
        delta_McFadden_R2_pp = delta_mcfadden_r2_pp
      ) %>%
      select(
        outcome,
        PGS,
        chosen_model,
        n_complete_cases,
        AUC_OOF_base,
        AUC_OOF_final,
        delta_AUC_OOF,
        Tjur_R2_covariates,
        Tjur_R2_selected,
        delta_Tjur_R2,
        delta_Tjur_R2_pp,
        McFadden_R2_covariates,
        McFadden_R2_selected,
        delta_McFadden_R2,
        delta_McFadden_R2_pp,
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
        pgs,
        sep = "__"
      )
    ]] <- tibble(
      outcome = outcome,
      PGS = pgs,
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
      delta_AUC_OOF = delta_auc_oof,
      Tjur_R2_covariates = tjur_r2_covariates,
      Tjur_R2_selected = tjur_r2_selected,
      delta_Tjur_R2 = delta_tjur_r2,
      delta_Tjur_R2_pp = delta_tjur_r2_pp,
      McFadden_R2_covariates = mcfadden_r2_covariates,
      McFadden_R2_selected = mcfadden_r2_selected,
      delta_McFadden_R2 = delta_mcfadden_r2,
      delta_McFadden_R2_pp = delta_mcfadden_r2_pp
    )

    model_performance_rows[[
      paste(
        outcome,
        pgs,
        sep = "__"
      )
    ]] <- tibble(
      outcome = outcome,
      PGS = pgs,
      chosen_model = chosen_model,
      selected_interactions = ifelse(
        length(final_interactions) == 0,
        "none",
        paste(final_interactions, collapse = ",")
      ),
      n_complete_cases = n_cc,
      AUC_OOF_selected = auc_oof_final,
      Tjur_R2_covariates = tjur_r2_covariates,
      Tjur_R2_selected = tjur_r2_selected,
      delta_Tjur_R2 = delta_tjur_r2,
      delta_Tjur_R2_pp = delta_tjur_r2_pp,
      McFadden_R2_covariates = mcfadden_r2_covariates,
      McFadden_R2_selected = mcfadden_r2_selected,
      delta_McFadden_R2 = delta_mcfadden_r2,
      delta_McFadden_R2_pp = delta_mcfadden_r2_pp
    )

    final_model_rows[[
      paste(
        outcome,
        pgs,
        sep = "__"
      )
    ]] <- final_tab

    final_models[[
      paste(
        outcome,
        pgs,
        sep = "__"
      )
    ]] <- list(
      fit = m_final,
      data = dat_cc,
      data_raw = dat_cc_raw,
      outcome = outcome,
      pgs = pgs,
      chosen_model = chosen_model,
      selected_interactions = final_interactions,
      site_terms = site_terms
    )

    oof_prediction_rows[[
      paste(
        outcome,
        pgs,
        "base",
        sep = "__"
      )
    ]] <- tibble(
      outcome = outcome,
      PGS = pgs,
      model = "base",
      y_true = dat_cc_raw[[outcome]],
      oof_pred = cv_base$oof_pred
    )

    oof_prediction_rows[[
      paste(
        outcome,
        pgs,
        "final",
        sep = "__"
      )
    ]] <- tibble(
      outcome = outcome,
      PGS = pgs,
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

model_performance_table <- bind_rows(
  model_performance_rows
)

oof_predictions_table <- bind_rows(
  oof_prediction_rows
)


# ============================================================
# 5. Retained interaction models
# ============================================================

selected_interaction_models <- model_selection_summary %>%
  filter(
    selected_interactions != "none"
  )


# ============================================================
# 6. Identify significant retained interaction terms
# ============================================================
# Only interaction terms present in the selected final models
# and showing Wald p < 0.05 are carried forward to post-hoc
# predicted-probability and simple-slope analyses.

significant_interactions <- final_OR_table %>%
  mutate(
    interaction = mapply(
      FUN = identify_interaction_type,
      term = term,
      pgs = PGS,
      MoreArgs = list(
        age_var = age_var,
        sex_var = sex_var,
        baseline_severity_var = baseline_severity_var
      ),
      USE.NAMES = FALSE
    )
  ) %>%
  filter(
    !is.na(interaction),
    p_value < 0.05
  ) %>%
  left_join(
    model_selection_summary %>%
      select(
        outcome,
        PGS,
        selected_interactions
      ),
    by = c(
      "outcome",
      "PGS"
    )
  ) %>%
  arrange(
    outcome,
    PGS,
    interaction
  )


# ============================================================
# 7. Audit retained interaction models
# ============================================================

interaction_model_audit_rows <- list()

for (
  i in seq_len(
    nrow(selected_interaction_models)
  )
) {

  outcome <- selected_interaction_models$outcome[i]
  pgs <- selected_interaction_models$PGS[i]

  key <- paste(
    outcome,
    pgs,
    sep = "__"
  )

  obj <- final_models[[key]]

  interaction_model_audit_rows[[
    length(interaction_model_audit_rows) + 1
  ]] <- tibble(
    outcome = outcome,
    PGS = pgs,
    selected_interactions = selected_interaction_models$selected_interactions[i],
    n_complete_cases = nrow(obj$data),
    formula = paste(
      deparse(
        formula(obj$fit)
      ),
      collapse = " "
    )
  )
}

interaction_model_audit <- bind_rows(
  interaction_model_audit_rows
)


# ============================================================
# 8. Predicted probabilities for significant retained
#    interactions
# ============================================================

predicted_probability_rows <- list()

for (
  i in seq_len(
    nrow(significant_interactions)
  )
) {

  outcome <- significant_interactions$outcome[i]
  pgs <- significant_interactions$PGS[i]
  interaction <- significant_interactions$interaction[i]

  key <- paste(
    outcome,
    pgs,
    sep = "__"
  )

  obj <- final_models[[key]]

  nd <- build_prediction_data(
    dat = obj$data,
    pgs = pgs,
    interaction = interaction,
    site_terms = obj$site_terms,
    pc_vars = pc_vars,
    age_var = age_var,
    sex_var = sex_var,
    baseline_severity_var = baseline_severity_var,
    pgs_values = prediction_pgs_values
  )

  nd <- add_prediction_intervals(
    fit = obj$fit,
    newdata = nd
  )

  nd <- nd %>%
    mutate(
      outcome = outcome,
      PGS = pgs,
      chosen_model = obj$chosen_model,
      interaction = interaction,
      interaction_term = significant_interactions$term[i],
      interaction_p_value = significant_interactions$p_value[i],
      interaction_OR = significant_interactions$OR[i]
    ) %>%
    select(
      outcome,
      PGS,
      chosen_model,
      interaction,
      interaction_term,
      interaction_p_value,
      interaction_OR,
      moderator_level,
      PGS_VALUE,
      predicted_prob,
      CI_low_prob,
      CI_high_prob,
      everything()
    )

  predicted_probability_rows[[
    paste(
      outcome,
      pgs,
      interaction,
      sep = "__"
    )
  ]] <- nd
}

predicted_probabilities <- bind_rows(
  predicted_probability_rows
)


# ============================================================
# 9. Simple-slope analyses and contrasts for significant
#    retained interactions
# ============================================================

simple_slope_rows <- list()
slope_contrast_rows <- list()

for (
  i in seq_len(
    nrow(significant_interactions)
  )
) {

  outcome <- significant_interactions$outcome[i]
  pgs <- significant_interactions$PGS[i]
  interaction <- significant_interactions$interaction[i]

  key <- paste(
    outcome,
    pgs,
    sep = "__"
  )

  obj <- final_models[[key]]

  fit <- obj$fit
  dat_cc <- obj$data
  dat_cc_raw <- obj$data_raw

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
      var = pgs
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
      pgs = pgs,
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
      pgs = pgs,
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
      raw_reference = dat_raw[[moderator]]
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
      var = pgs,
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
      pgs = pgs,
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
      pgs = pgs,
      outcome = outcome,
      moderator = moderator,
      model_n = nrow(dat_cc),
      contrast_labels = contrast_labels
    )
  }
}

simple_pgs_slopes <- bind_rows(
  simple_slope_rows
)

pgs_slope_contrasts <- bind_rows(
  slope_contrast_rows
)


# ============================================================
# 10. Save outputs
# ============================================================


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
  model_performance_table,
  file.path(
    output_dir,
    "model_performance_table.csv"
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
  significant_interactions,
  file.path(
    output_dir,
    "significant_interactions.csv"
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
  simple_pgs_slopes,
  file.path(
    output_dir,
    "simple_pgs_slopes.csv"
  )
)

write_csv(
  pgs_slope_contrasts,
  file.path(
    output_dir,
    "pgs_slope_contrasts.csv"
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
    model_performance = model_performance_table,
    OOF_predictions = oof_predictions_table,
    selected_interactions = selected_interaction_models,
    significant_interactions = significant_interactions,
    predicted_probabilities = predicted_probabilities,
    simple_slopes = simple_pgs_slopes,
    slope_contrasts = pgs_slope_contrasts,
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
