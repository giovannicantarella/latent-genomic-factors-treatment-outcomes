# ============================================================
# Project: Latent genomic factors and treatment outcomes
# Script: Preparation of a cohort dataset with PRS and covariates
#
# Purpose:
# - read a cohort-level dataset
# - read and merge polygenic score files
# - merge external covariates
# - harmonise duplicated columns created during merging
# - convert numeric-like character variables where needed
# - derive analysis-ready outcomes when required
# - save the final analysis-ready dataset
#
# Notes:
# - this script is written as a generic template and should be
#   adapted to the structure of the cohort being processed
# - file paths, separators, subject identifiers, and outcome
#   definitions may differ across datasets
# - subject matching is assumed to rely on one identifier shared
#   across the cohort dataset, PRS files, and covariate files,
#   after harmonisation if needed
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# ============================================================
# 1. Input files and directories
# ============================================================
# Replace these placeholders with the relevant paths for the
# cohort of interest.
# ============================================================

prs_dir         <- "/path/to/directory/containing/sscore/files"
cohort_file     <- "/path/to/cohort_dataset.csv"
covariates_file <- "/path/to/external_covariates.txt"
output_file     <- "/path/to/output/cohort_prepared_with_PRS.csv"

# Update this separator if needed
cohort_sep <- ";"

# ============================================================
# 2. Helper functions
# ============================================================

comma_to_numeric <- function(x) {
  if (!is.character(x)) return(x)
  x <- trimws(x)
  x[x == ""] <- NA
  suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
}

is_num_like <- function(v) {
  is.character(v) &&
    all(is.na(v) | grepl("^\\s*[+-]?\\d+(?:[\\.,]\\d+)?\\s*$", v))
}

safe_z <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

canon_id <- function(x) {
  toupper(trimws(as.character(x)))
}

coalesce_xy <- function(df, prefix_regex) {
  xy <- grep(paste0("^", prefix_regex, "\\.(x|y)$"), names(df), value = TRUE)
  if (length(xy) == 0) return(df)
  
  base <- unique(str_remove(xy, "\\.(x|y)$"))
  for (b in base) {
    xcol <- paste0(b, ".x")
    ycol <- paste0(b, ".y")
    
    if (xcol %in% names(df) && ycol %in% names(df)) {
      df[[b]] <- dplyr::coalesce(df[[xcol]], df[[ycol]])
    } else if (xcol %in% names(df)) {
      df[[b]] <- df[[xcol]]
    } else if (ycol %in% names(df)) {
      df[[b]] <- df[[ycol]]
    }
  }
  
  rm_cols <- intersect(names(df), c(paste0(base, ".x"), paste0(base, ".y")))
  df %>% select(-any_of(rm_cols))
}

read_and_rename_sscore <- function(path) {
  filename <- basename(path)
  
  if (grepl("^PRScs_.+?_PRScs_scores\\.sscore$", filename)) {
    trait <- sub("^PRScs_(.+?)_PRScs_scores\\.sscore$", "\\1", filename)
  } else if (grepl("^PRScs_scores_.+?\\.sscore$", filename)) {
    trait <- sub("^PRScs_scores_(.+?)\\.sscore$", "\\1", filename)
  } else {
    stop("Unrecognised .sscore filename format: ", filename)
  }
  
  df <- read_tsv(path, show_col_types = FALSE)
  
  required_cols <- c("IID", "SCORE1_AVG")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", filename, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    select(IID, SCORE1_AVG) %>%
    rename(!!paste0("PRS_", trait) := SCORE1_AVG)
}

is_binary_int <- function(v) {
  is.integer(v) && all(is.na(v) | v %in% c(0L, 1L))
}

# ============================================================
# 3. Read cohort dataset
# ============================================================
# Replace GENETIC_ID with the identifier used in the cohort-
# level dataset.
# ============================================================

cohort_df <- read.csv(cohort_file, sep = cohort_sep, stringsAsFactors = FALSE)

if (!"GENETIC_ID" %in% names(cohort_df)) {
  stop("The cohort dataset must contain a 'GENETIC_ID' column.")
}

cohort_df$GENETIC_ID <- canon_id(cohort_df$GENETIC_ID)

# ============================================================
# 4. Read and merge PRS files
# ============================================================
# PRS files are expected to be stored as .sscore files.
# Replace the merging key below if the cohort uses a different
# identifier mapping across the cohort dataset, PRS files,
# and covariate files.
# ============================================================

sscore_paths <- list.files(
  path = prs_dir,
  pattern = "\\.sscore$",
  full.names = TRUE
)

if (length(sscore_paths) == 0) {
  stop("No .sscore files were found in the specified directory.")
}

prs_list <- lapply(sscore_paths, read_and_rename_sscore)
prs_wide <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), prs_list)

prs_wide$IID <- canon_id(prs_wide$IID)

if (anyDuplicated(prs_wide$IID)) {
  stop("Duplicated IID values were found in the merged PRS dataset.")
}

cohort_df <- merge(cohort_df, prs_wide, by.x = "GENETIC_ID", by.y = "IID", all.x = TRUE)

# ============================================================
# 5. Add PRS z-scores
# ============================================================

prs_cols <- grep("^PRS_", names(cohort_df), value = TRUE)

for (col in prs_cols) {
  cohort_df[[paste0(col, "_z")]] <- safe_z(cohort_df[[col]])
}

# Optional matching checks
ids_prs     <- prs_wide$IID
ids_genetic <- cohort_df$GENETIC_ID

cat("PRS IDs without match in cohort dataset:", sum(!ids_prs %in% ids_genetic, na.rm = TRUE), "\n")
cat("Cohort IDs without match in PRS:", sum(!ids_genetic %in% ids_prs, na.rm = TRUE), "\n")

cohort_df <- cohort_df[!duplicated(cohort_df$GENETIC_ID), , drop = FALSE]

# ============================================================
# 6. Read and merge external covariates
# ============================================================
# Replace the expected identifier name if the external covariate
# file uses a different column.
# ============================================================

covariates_df <- read.table(
  covariates_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

if ("#IID" %in% names(covariates_df)) {
  names(covariates_df)[names(covariates_df) == "#IID"] <- "IID"
}

if (!"IID" %in% names(covariates_df)) {
  stop("The covariate file must contain an 'IID' column (or '#IID').")
}

covariates_df$IID <- canon_id(covariates_df$IID)

dup_n <- sum(duplicated(covariates_df$IID))
if (dup_n > 0) {
  message(
    "Duplicated IID values found in the covariate file (keeping the first occurrence): ",
    dup_n
  )
  covariates_df <- covariates_df[!duplicated(covariates_df$IID), ]
}

cohort_df <- merge(cohort_df, covariates_df, by.x = "GENETIC_ID", by.y = "IID", all.x = TRUE)

cat("Rows after covariate merge:", nrow(cohort_df), "\n")

# ============================================================
# 7. Harmonise duplicated columns after merging
# ============================================================
# Adjust these patterns depending on which duplicated variables
# are expected after merging in the cohort of interest.
# ============================================================

cohort_df <- coalesce_xy(cohort_df, "PC\\d+")
cohort_df <- coalesce_xy(cohort_df, "site\\d+")
cohort_df <- coalesce_xy(cohort_df, "age")
cohort_df <- coalesce_xy(cohort_df, "Age")
cohort_df <- coalesce_xy(cohort_df, "sex")
cohort_df <- coalesce_xy(cohort_df, "Sex")
cohort_df <- coalesce_xy(cohort_df, "outcome")
cohort_df <- coalesce_xy(cohort_df, "remission")
cohort_df <- coalesce_xy(cohort_df, "response")
cohort_df <- coalesce_xy(cohort_df, "resistance")

if (!"age" %in% names(cohort_df) && "Age" %in% names(cohort_df)) {
  cohort_df$age <- cohort_df$Age
}

if (!"sex" %in% names(cohort_df) && "Sex" %in% names(cohort_df)) {
  cohort_df$sex <- cohort_df$Sex
}

# ============================================================
# 8. Convert numeric-like columns and set variable types
# ============================================================

num_like_cols <- names(cohort_df)[vapply(cohort_df, is_num_like, logical(1))]
if (length(num_like_cols) > 0) {
  cohort_df <- cohort_df %>%
    mutate(across(all_of(num_like_cols), comma_to_numeric))
}

pc_cols <- grep("^PC[0-9]+$", names(cohort_df), value = TRUE)
if (length(pc_cols) > 0) {
  cohort_df <- cohort_df %>%
    mutate(across(all_of(pc_cols), as.numeric))
}

site_cols <- intersect(paste0("site", 1:20), names(cohort_df))
if (length(site_cols) > 0) {
  cohort_df <- cohort_df %>%
    mutate(across(all_of(site_cols), ~ factor(., levels = c(0, 1))))
}

bin_cols <- setdiff(
  names(cohort_df)[vapply(cohort_df, is_binary_int, logical(1))],
  site_cols
)

if (length(bin_cols) > 0) {
  cohort_df <- cohort_df %>%
    mutate(across(all_of(bin_cols), ~ factor(., levels = c(0, 1))))
}

# ============================================================
# 9. Derive binary outcomes where needed
# ============================================================
# This section is intentionally generic.
# Replace SOURCE_OUTCOME and the coding logic below with the
# definitions used in the cohort being processed.
# ============================================================

if ("SOURCE_OUTCOME" %in% names(cohort_df)) {
  cohort_df$Response <- ifelse(cohort_df$SOURCE_OUTCOME == 1, 1, 0)
  cohort_df$Resistance <- ifelse(cohort_df$SOURCE_OUTCOME == 3, 1, 0)
  
  cat("\nResponse (derived outcome):\n")
  print(table(cohort_df$Response, useNA = "ifany"))
  
  cat("\nResistance (derived outcome):\n")
  print(table(cohort_df$Resistance, useNA = "ifany"))
}

# ============================================================
# 10. Check selected column classes
# ============================================================

prs_cols_all <- grep("^PRS_", names(cohort_df), value = TRUE)
prs_z <- prs_cols_all[str_detect(prs_cols_all, "_z$")]
prs_check <- if (length(prs_z) > 0) prs_z else prs_cols_all

check_cols <- c("CURRENT_SEVERITY", prs_check)
check_cols <- intersect(check_cols, names(cohort_df))

if (length(check_cols) > 0) {
  classes <- sapply(cohort_df[check_cols], function(x) paste(class(x), collapse = "/"))
  print(classes)
}

# ============================================================
# 11. Save final dataset
# ============================================================

write.csv(cohort_df, output_file, row.names = FALSE)

cat("\nSaved prepared dataset to:\n", output_file, "\n", sep = "")
cat("Final number of rows:", nrow(cohort_df), "\n")
cat("Final number of columns:", ncol(cohort_df), "\n")