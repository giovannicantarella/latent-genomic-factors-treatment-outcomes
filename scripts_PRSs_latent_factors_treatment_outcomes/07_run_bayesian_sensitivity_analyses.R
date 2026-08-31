# ============================================================
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Bayesian sensitivity analysis
#
# Purpose:
# - read the models selected in the frequentist analysis
# - refit the same models using Bayesian logistic regression
# - apply weakly informative priors to regularise coefficients
# - estimate posterior odds ratios and 95% credible intervals
# - calculate posterior predicted probabilities for retained
#   interaction models
# - save Bayesian sensitivity-analysis outputs
#
# Notes:
# - model selection is not repeated in this script
# - each Bayesian model reproduces the corresponding model
#   selected in the frequentist analysis
# - regression coefficients use Normal(0, 0.5) priors on the
#   log-odds scale, with autoscaling disabled
# - the intercept uses a weakly informative Student-t prior
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(rstanarm)
library(posterior)
library(readr)
library(writexl)


# ============================================================
# 1. User-defined settings
# ============================================================

data_file <- "analysis_ready_dataset.csv"

# Output from 05_run_frequentist_models.R
selection_file <- file.path(
  "frequentist_outputs",
  "model_selection_summary.csv"
)

output_dir <- "bayesian_sensitivity_outputs"

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

site_terms_by_outcome <- list(
  OUTCOME_1 = c("SITE_1", "SITE_2"),
  OUTCOME_2 = c("SITE_1", "SITE_2"),
  OUTCOME_3 = c("SITE_1", "SITE_2"),
  OUTCOME_4 = c("SITE_1", "SITE_2"),
  OUTCOME_5 = c("SITE_1", "SITE_2")
)

# MCMC settings.
# These values should be chosen so that posterior sampling is stable;
# convergence diagnostics should be checked before interpreting
# posterior estimates.
stan_chains <- 4
stan_iter <- 2000
stan_seed <- 123
stan_refresh <- 100

options(
  mc.cores = min(
    stan_chains,
    parallel::detectCores()
  )
)

# Predicted probabilities use the same PGS contrast as the
# frequentist post-hoc analyses. Because the PGS is Gelman-scaled,
# -2 and +2 SD on the original scale correspond to -1 and +1 on
# the model scale.
prediction_prs_values <- c(-1, 1)
moderator_probs <- c(0.25, 0.50, 0.75)


# ============================================================
# 2. Helper functions
# ============================================================

# Apply the same Gelman scaling used in the frequentist analysis so
# that the Bayesian models reproduce the same model specification.
gelman_scale <- function(x) {
  if (all(is.na(x))) return(x)

  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(x - m)
  }

  (x - m) / (2 * s)
}


# Verify that binary outcomes use 0/1 coding before model fitting.
# The script stops if an outcome uses any other coding.
as_binary_01 <- function(x, name = "outcome") {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    x <- suppressWarnings(
      as.numeric(x)
    )
  }

  values <- sort(
    unique(
      x[!is.na(x)]
    )
  )

  if (!all(values %in% c(0, 1))) {
    stop(
      name,
      " must be coded 0/1. Observed values: ",
      paste(
        values,
        collapse = ", "
      )
    )
  }

  as.integer(x)
}


# Reconstruct the model formula selected in the frequentist analysis.
# No new interaction search or model comparison is performed here.
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
      paste0(
        prs,
        ":",
        baseline_severity_var
      )
    )
  }

  if ("sex" %in% interactions) {
    rhs_terms <- c(
      rhs_terms,
      paste0(
        prs,
        ":",
        sex_var
      )
    )
  }

  if ("age" %in% interactions) {
    rhs_terms <- c(
      rhs_terms,
      paste0(
        prs,
        ":",
        age_var
      )
    )
  }

  rhs_string <- paste(
    c(
      rhs_terms,
      site_terms
    ),
    collapse = " + "
  )

  as.formula(
    paste(
      outcome,
      "~",
      rhs_string
    )
  )
}


# Recreate the complete-case dataset used for a given PGS-outcome
# model using the same covariates and site terms as the frequentist
# analysis.
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
    filter(
      !is.na(
        .data[[outcome]]
      )
    ) %>%
    mutate(
      "{sex_var}" := factor(
        .data[[sex_var]]
      )
    ) %>%
    select(
      all_of(
        needed_vars
      )
    ) %>%
    drop_na()
}


# Convert the interaction labels written by the frequentist script
# back into a character vector used to rebuild the selected formula.
interaction_from_selection <- function(x) {
  if (
    is.na(x) ||
      x == "none"
  ) {
    return(
      character(0)
    )
  }

  unlist(
    strsplit(
      x,
      "\\s*,\\s*"
    )
  )
}


# Summarise posterior coefficient draws.
# Posterior means and 95% equal-tail credible intervals are reported
# on both the log-odds and odds-ratio scales.
tidy_bayes_or <- function(fit) {

  coefficient_names <- names(
    coef(fit)
  )

  draws <- as.matrix(
    fit
  )

  draws <- draws[
    ,
    intersect(
      colnames(draws),
      coefficient_names
    ),
    drop = FALSE
  ]

  estimate <- apply(
    draws,
    2,
    mean
  )

  posterior_sd <- apply(
    draws,
    2,
    sd
  )

  q025 <- apply(
    draws,
    2,
    quantile,
    probs = 0.025,
    names = FALSE
  )

  q975 <- apply(
    draws,
    2,
    quantile,
    probs = 0.975,
    names = FALSE
  )

  tibble(
    term = names(estimate),
    estimate = as.numeric(estimate),
    posterior_sd = as.numeric(posterior_sd),
    CI_low_beta = as.numeric(q025),
    CI_high_beta = as.numeric(q975),
    OR = exp(estimate),
    CI_low = exp(q025),
    CI_high = exp(q975)
  )
}


# Use the same covariate reference profile as in the frequentist
# probability calculations so that the two sets of predictions are
# directly comparable.
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


# Build the design matrix for new prediction data using the terms and
# contrasts stored in the fitted Bayesian model.
model_matrix_from_fit <- function(
  fit,
  newdata
) {

  model_terms <- delete.response(
    terms(fit)
  )

  model.matrix(
    model_terms,
    data = newdata,
    contrasts.arg = fit$contrasts
  )
}


# Convert posterior linear predictors to probabilities for each draw,
# then summarise the posterior probability distribution using its mean
# and 95% credible interval.
posterior_prediction_summary <- function(
  fit,
  newdata
) {

  coefficient_names <- names(
    coef(fit)
  )

  draws <- as.matrix(
    fit
  )

  draws <- draws[
    ,
    intersect(
      colnames(draws),
      coefficient_names
    ),
    drop = FALSE
  ]

  X <- model_matrix_from_fit(
    fit,
    newdata
  )

  common_terms <- intersect(
    colnames(X),
    colnames(draws)
  )

  X <- X[
    ,
    common_terms,
    drop = FALSE
  ]

  draws <- draws[
    ,
    common_terms,
    drop = FALSE
  ]

  eta <- draws %*% t(X)
  probability <- plogis(eta)

  newdata %>%
    mutate(
      predicted_prob = apply(
        probability,
        2,
        mean,
        na.rm = TRUE
      ),
      CrI_low_prob = apply(
        probability,
        2,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      CrI_high_prob = apply(
        probability,
        2,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      )
    )
}


# Reproduce the same post-hoc prediction grid used in the frequentist
# analysis: P25/P50/P75 for continuous moderators, observed levels
# for sex, and -2/+2 raw SD for the PGS.
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
      SEX_LEVEL = levels(
        dat[[sex_var]]
      ),
      stringsAsFactors = FALSE
    )

    nd[[prs]] <- nd$PRS_VALUE

    nd[[age_var]] <- ref[[age_var]]
    nd[[baseline_severity_var]] <- ref[[baseline_severity_var]]

    nd[[sex_var]] <- factor(
      nd$SEX_LEVEL,
      levels = levels(
        dat[[sex_var]]
      )
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
      levels = levels(
        dat[[sex_var]]
      )
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
      levels = levels(
        dat[[sex_var]]
      )
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


# ============================================================
# 3. Read and prepare data
# ============================================================

# The frequentist model-selection table is treated as the source of
# the final model specification. Bayesian analysis is therefore a
# sensitivity re-estimation step rather than a second selection stage.
dat_raw <- read.csv(
  data_file,
  stringsAsFactors = FALSE
)

model_selection <- read_csv(
  selection_file,
  show_col_types = FALSE
)

dat_raw <- dat_raw %>%
  filter(
    is.na(
      .data[[age_var]]
    ) |
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

dat_scaled <- dat_raw %>%
  mutate(
    across(
      all_of(
        existing_continuous_vars
      ),
      gelman_scale
    )
  )

model_map <- model_selection %>%
  select(
    outcome,
    PGS,
    chosen_model,
    selected_interactions,
    n_complete_cases,
    AUC_OOF_base,
    AUC_OOF_final,
    delta_AUC_OOF
  ) %>%
  distinct()

duplicate_check <- model_map %>%
  count(
    outcome,
    PGS
  )

if (
  any(
    duplicate_check$n > 1
  )
) {
  stop(
    "The frequentist model-selection file contains duplicate outcome x PGS rows."
  )
}


# ============================================================
# 4. Fit Bayesian versions of selected models
# ============================================================

# Each selected frequentist model is refitted with stan_glm().
#
# Normal(0, 0.5) priors are assigned to regression coefficients on
# the log-odds scale. These priors favour smaller coefficient values
# and pull weakly supported estimates toward zero without fixing them
# at the null.
#
# A separate weakly informative Student-t prior is used for the
# intercept. autoscale = FALSE keeps the stated prior scales fixed.
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

bayesian_summary_rows <- list()
bayesian_model_rows <- list()
bayesian_models <- list()

for (outcome in outcomes) {

  site_terms <- site_terms_by_outcome[[outcome]]

  if (is.null(site_terms)) {
    stop(
      "No site-variable specification found for ",
      outcome
    )
  }

  for (prs in prs_vars) {

    freq_row <- model_map %>%
      filter(
        .data$outcome == outcome,
        .data$PGS == prs
      )

    if (nrow(freq_row) == 0) {
      warning(
        "No frequentist selected model found for ",
        outcome,
        " and ",
        prs
      )
      next
    }

    interactions <- interaction_from_selection(
      freq_row$selected_interactions[1]
    )

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

    n_cc <- nrow(
      dat_cc
    )

    if (n_cc == 0) {
      warning(
        "No complete-case data for ",
        outcome,
        " and ",
        prs
      )
      next
    }

    final_formula <- build_formula_from_interactions(
      outcome = outcome,
      prs = prs,
      interactions = interactions,
      site_terms = site_terms,
      pc_vars = pc_vars,
      age_var = age_var,
      sex_var = sex_var,
      baseline_severity_var = baseline_severity_var
    )

    # Refit the exact frequentist-selected formula. This preserves
    # the same outcome definition, complete-case sample, covariates,
    # site adjustment, and retained interaction terms.
    fit <- stan_glm(
      final_formula,
      data = dat_cc,
      family = binomial(
        link = "logit"
      ),
      prior = normal(
        location = 0,
        scale = 0.5,
        autoscale = FALSE
      ),
      prior_intercept = student_t(
        df = 7,
        location = 0,
        scale = 2.5,
        autoscale = FALSE
      ),
      chains = stan_chains,
      iter = stan_iter,
      seed = stan_seed,
      refresh = stan_refresh
    )

    bayes_tab <- tidy_bayes_or(
      fit
    ) %>%
      mutate(
        outcome = outcome,
        PGS = prs,
        chosen_model = freq_row$chosen_model[1],
        selected_interactions = freq_row$selected_interactions[1],
        n_complete_cases = n_cc,
        n_complete_cases_freq = freq_row$n_complete_cases[1],
        AUC_OOF_base_freq = freq_row$AUC_OOF_base[1],
        AUC_OOF_final_freq = freq_row$AUC_OOF_final[1],
        delta_AUC_OOF_freq = freq_row$delta_AUC_OOF[1]
      ) %>%
      select(
        outcome,
        PGS,
        chosen_model,
        selected_interactions,
        n_complete_cases,
        n_complete_cases_freq,
        AUC_OOF_base_freq,
        AUC_OOF_final_freq,
        delta_AUC_OOF_freq,
        term,
        estimate,
        posterior_sd,
        OR,
        CI_low,
        CI_high,
        CI_low_beta,
        CI_high_beta
      )

    bayesian_model_rows[[
      paste(
        outcome,
        prs,
        sep = "__"
      )
    ]] <- bayes_tab

    bayesian_summary_rows[[
      paste(
        outcome,
        prs,
        sep = "__"
      )
    ]] <- tibble(
      outcome = outcome,
      PGS = prs,
      chosen_model = freq_row$chosen_model[1],
      selected_interactions = freq_row$selected_interactions[1],
      n_complete_cases = n_cc,
      n_complete_cases_freq = freq_row$n_complete_cases[1],
      AUC_OOF_base_freq = freq_row$AUC_OOF_base[1],
      AUC_OOF_final_freq = freq_row$AUC_OOF_final[1],
      delta_AUC_OOF_freq = freq_row$delta_AUC_OOF[1],
      site_terms = paste(
        site_terms,
        collapse = " + "
      )
    )

    bayesian_models[[
      paste(
        outcome,
        prs,
        sep = "__"
      )
    ]] <- list(
      fit = fit,
      data = dat_cc,
      outcome = outcome,
      prs = prs,
      chosen_model = freq_row$chosen_model[1],
      selected_interactions = interactions,
      site_terms = site_terms
    )
  }
}

bayesian_model_summary <- bind_rows(
  bayesian_summary_rows
)

bayesian_OR_table <- bind_rows(
  bayesian_model_rows
)


# ============================================================
# 5. Posterior predicted probabilities for retained interactions
# ============================================================

# Posterior probabilities are calculated only for models containing
# retained interactions. The prediction grid matches the frequentist
# analysis, allowing the fitted interaction patterns and uncertainty
# intervals to be compared across the two approaches.
bayesian_prediction_rows <- list()

interaction_models <- bayesian_model_summary %>%
  filter(
    selected_interactions != "none"
  )

for (
  i in seq_len(
    nrow(interaction_models)
  )
) {

  outcome <- interaction_models$outcome[i]
  prs <- interaction_models$PGS[i]

  interactions <- interaction_from_selection(
    interaction_models$selected_interactions[i]
  )

  key <- paste(
    outcome,
    prs,
    sep = "__"
  )

  obj <- bayesian_models[[key]]

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

    nd <- posterior_prediction_summary(
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
        CrI_low_prob,
        CrI_high_prob,
        everything()
      )

    bayesian_prediction_rows[[
      paste(
        outcome,
        prs,
        interaction,
        sep = "__"
      )
    ]] <- nd
  }
}

bayesian_predicted_probabilities <- bind_rows(
  bayesian_prediction_rows
)


# ============================================================
# 6. Save outputs
# ============================================================

# Save model-level information, posterior coefficient estimates, and
# model-based probabilities as separate CSV files and in one workbook.
write_csv(
  bayesian_model_summary,
  file.path(
    output_dir,
    "bayesian_model_summary.csv"
  )
)

write_csv(
  bayesian_OR_table,
  file.path(
    output_dir,
    "bayesian_OR_table.csv"
  )
)

write_csv(
  bayesian_predicted_probabilities,
  file.path(
    output_dir,
    "bayesian_predicted_probabilities.csv"
  )
)

write_xlsx(
  list(
    model_summary = bayesian_model_summary,
    posterior_ORs = bayesian_OR_table,
    predicted_probabilities = bayesian_predicted_probabilities
  ),
  path = file.path(
    output_dir,
    "bayesian_sensitivity_outputs.xlsx"
  )
)

cat(
  "\nBayesian sensitivity-analysis outputs saved in:\n",
  output_dir,
  "\n"
)
