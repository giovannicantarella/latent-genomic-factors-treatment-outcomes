# ============================================================
# Prepare GSRD dataset with PRS and covariates
# ============================================================
# This script:
# - reads the GSRD clinical dataset
# - reads all PRS-CS .sscore files from a specified directory
# - merges PRS scores into the clinical dataset
# - computes z-scores for PRS variables
# - merges external covariates
# - harmonises duplicated columns created during merging
# - converts numeric-like character variables to numeric
# - creates binary outcomes where needed
# - saves the final analysis-ready dataset
#
# Notes:
# - Update all file paths before running
# - Subject matching is based on GenID (clinical data) and IID
#   (PRS/covariates)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# ============================================================
# 1. Input files and directories
# ============================================================

# Insert your paths here
prs_dir         <- "/path/to/directory/containing/sscore/files"
gsrd_file       <- "/path/to/GSRD.csv"
covariates_file <- "/path/to/gsrd_cov.txt"
output_file     <- "/path/to/output/GSRD_prepared_with_PRS.csv"

# Update this separator if needed
gsrd_sep <- ";"

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
    stop("Missing required columns in ", filename, ": ",
         paste(missing_cols, collapse = ", "))
  }
  
  df %>%
    select(IID, SCORE1_AVG) %>%
    rename(!!paste0("PRS_", trait) := SCORE1_AVG)
}

# ============================================================
# 3. Read GSRD clinical dataset
# ============================================================

GSRD <- read.csv(gsrd_file, sep = gsrd_sep, stringsAsFactors = FALSE)

if (!"GenID" %in% names(GSRD)) {
  stop("The GSRD dataset must contain a 'GenID' column.")
}

GSRD$GenID <- toupper(trimws(as.character(GSRD$GenID)))

# ============================================================
# 4. Read and merge PRS files
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

if (anyDuplicated(prs_wide$IID)) {
  stop("Duplicated IID values were found in the merged PRS dataset.")
}

GSRD <- merge(GSRD, prs_wide, by.x = "GenID", by.y = "IID", all.x = TRUE)

# ============================================================
# 5. Add PRS z-scores
# ============================================================

prs_cols <- grep("^PRS_", names(GSRD), value = TRUE)

for (col in prs_cols) {
  GSRD[[paste0(col, "_z")]] <- safe_z(GSRD[[col]])
}

# Optional matching checks
ids_prs  <- prs_wide$IID
ids_gsrd <- GSRD$GenID

cat("PRS IDs without match in GSRD:", sum(!ids_prs %in% ids_gsrd, na.rm = TRUE), "\n")
cat("GSRD IDs without match in PRS:", sum(!ids_gsrd %in% ids_prs, na.rm = TRUE), "\n")

# Remove duplicated GenID if present
GSRD <- GSRD[!duplicated(GSRD$GenID), , drop = FALSE]

# ============================================================
# 6. Read and merge external covariates
# ============================================================

gsrd_cov <- read.table(
  covariates_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

if ("#IID" %in% names(gsrd_cov)) {
  names(gsrd_cov)[names(gsrd_cov) == "#IID"] <- "IID"
}

if (!"IID" %in% names(gsrd_cov)) {
  stop("The covariate file must contain an 'IID' column (or '#IID').")
}

gsrd_cov$IID <- toupper(trimws(as.character(gsrd_cov$IID)))

dup_n <- sum(duplicated(gsrd_cov$IID))
if (dup_n > 0) {
  message(
    "Duplicated IID values found in the covariate file (keeping the first occurrence): ",
    dup_n
  )
  gsrd_cov <- gsrd_cov[!duplicated(gsrd_cov$IID), ]
}

GSRD <- merge(GSRD, gsrd_cov, by.x = "GenID", by.y = "IID", all.x = TRUE)

cat("Rows after covariate merge:", nrow(GSRD), "\n")

# ============================================================
# 7. Harmonise duplicated columns after merging
# ============================================================

GSRD <- coalesce_xy(GSRD, "PC\\d+")
GSRD <- coalesce_xy(GSRD, "site\\d+")
GSRD <- coalesce_xy(GSRD, "age")
GSRD <- coalesce_xy(GSRD, "Age")
GSRD <- coalesce_xy(GSRD, "sex")
GSRD <- coalesce_xy(GSRD, "Sex")
GSRD <- coalesce_xy(GSRD, "MADRS\\.Retrospective")
GSRD <- coalesce_xy(GSRD, "MADRS\\.Current")
GSRD <- coalesce_xy(GSRD, "MADRS\\.Outcome")
GSRD <- coalesce_xy(GSRD, "response")
GSRD <- coalesce_xy(GSRD, "remission")
GSRD <- coalesce_xy(GSRD, "resistance")
GSRD <- coalesce_xy(GSRD, "augmentation_antipsych")
GSRD <- coalesce_xy(GSRD, "augmentation_moodstab")

if (!"age" %in% names(GSRD) && "Age" %in% names(GSRD)) {
  GSRD$age <- GSRD$Age
}

if (!"sex" %in% names(GSRD) && "Sex" %in% names(GSRD)) {
  GSRD$sex <- GSRD$Sex
}

# ============================================================
# 8. Convert numeric-like columns and set variable types
# ============================================================

num_like_cols <- names(GSRD)[vapply(GSRD, is_num_like, logical(1))]
if (length(num_like_cols) > 0) {
  GSRD <- GSRD %>%
    mutate(across(all_of(num_like_cols), comma_to_numeric))
}

pc_cols <- grep("^PC[0-9]+$", names(GSRD), value = TRUE)
if (length(pc_cols) > 0) {
  GSRD <- GSRD %>%
    mutate(across(all_of(pc_cols), as.numeric))
}

site_cols <- intersect(paste0("site", 1:9), names(GSRD))
if (length(site_cols) > 0) {
  GSRD <- GSRD %>%
    mutate(across(all_of(site_cols), ~ factor(., levels = c(0, 1))))
}

is_binary_int <- function(v) {
  is.integer(v) && all(is.na(v) | v %in% c(0L, 1L))
}

bin_cols <- setdiff(
  names(GSRD)[vapply(GSRD, is_binary_int, logical(1))],
  site_cols
)

if (length(bin_cols) > 0) {
  GSRD <- GSRD %>%
    mutate(across(all_of(bin_cols), ~ factor(., levels = c(0, 1))))
}

# ============================================================
# 9. Create binary outcomes from MADRS.Outcome, if available
# ============================================================

if ("MADRS.Outcome" %in% names(GSRD)) {
  GSRD$Response   <- ifelse(GSRD$MADRS.Outcome == 1, 0, 1)
  GSRD$Resistance <- ifelse(GSRD$MADRS.Outcome == 3, 1, 0)
  
  cat("\nResponse (derived from MADRS.Outcome):\n")
  print(table(GSRD$Response, useNA = "ifany"))
  
  cat("\nResistance (derived from MADRS.Outcome):\n")
  print(table(GSRD$Resistance, useNA = "ifany"))
}

# ============================================================
# 10. Check selected column classes
# ============================================================

prs_cols_all <- grep("^PRS_", names(GSRD), value = TRUE)
prs_z <- prs_cols_all[str_detect(prs_cols_all, "_z$")]
prs_check <- if (length(prs_z) > 0) prs_z else prs_cols_all

check_cols <- c("MADRS.Current", prs_check)
check_cols <- intersect(check_cols, names(GSRD))

if (length(check_cols) > 0) {
  classes <- sapply(GSRD[check_cols], function(x) paste(class(x), collapse = "/"))
  print(classes)
}

# ============================================================
# 11. Save final dataset
# ============================================================

write.csv(GSRD, output_file, row.names = FALSE)

cat("\nSaved prepared dataset to:\n", output_file, "\n", sep = "")
cat("Final number of rows:", nrow(GSRD), "\n")
cat("Final number of columns:", ncol(GSRD), "\n")