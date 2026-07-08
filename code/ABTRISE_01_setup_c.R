# =============================================================================
# ABT-RISE: Site-Level Analysis Script 1 of 4
# SETUP -- Configuration, Data Load, Preparation, Shared Diagnostics
#
# SCRIPTS IN THIS SERIES:
#   ABTRISE_01_setup.R          <- YOU ARE HERE
#   ABTRISE_02_criterion.R      source("ABTRISE_01_setup.R") automatically
#   ABTRISE_345_outcomes.R      source("ABTRISE_01_setup.R") automatically
#   ABTRISE_06_benchmarking.R   source("ABTRISE_01_setup.R") automatically
#
# HOW TO RUN:
#   - Run this script once directly to check data quality and review
#     diagnostics before running any analysis script.
#   - Each analysis script (02, 345, 06) sources this file automatically
#     at startup -- you do not need to run this manually before them.
#   - All four scripts are independent and can be run in any order.
#     Each re-runs setup fresh to ensure a clean environment.
#
# WHAT THIS SCRIPT PRODUCES:
#   outputs/diagnostics/  cohort_summary diagnostics, missingness,
#                         exclusion waterfall, impossible values,
#                         delivery rate tables and figures,
#                         table 1, sedation distribution
#
# DATA INPUTS:
#   File 1: file1_person_period.parquet   -- one row per patient per vent-day
#   File 2: file2_hospitalization_level.parquet -- one row per hospitalization
#
# COORDINATING CENTER: Rush
# =============================================================================

# =============================================================================
# SECTION 0: SITE CONFIGURATION
# =============================================================================
# *** SITES: ONLY EDIT THIS SECTION ***
# Everything below Section 0 runs automatically.

# Site config (site_id, data_dir, out_dir) comes from clif_config.json via
# ABTRISE_config.R. run_all.R sources it; if this script is sourced directly
# from RStudio, load the config here as a fallback.
if (!exists("site_id") || !exists("data_dir") || !exists("out_dir")) {
  suppressPackageStartupMessages(library(here))
  source(here::here("code", "ABTRISE_config.R"))
}

# =============================================================================
# SECTION 1: LIBRARIES, HELPERS, OUTPUT FOLDERS
# =============================================================================

# ensure_packages() is defined in ABTRISE_config.R, which is always sourced
# before this block (by run_all.R, or by the fallback at lines 41-44 above).
# Install any missing packages before the fail-fast library() calls below.
# `splines` is base R (ships with R) and is intentionally omitted.
ensure_packages(c(
  "here", "arrow", "dplyr", "tidyr", "stringr", "lme4", "glmmTMB",
  "survival", "tidycmprsk", "broom.mixed", "broom", "readr", "ggplot2",
  "patchwork", "purrr", "scales", "epiR", "blandr", "forcats", "flextable"
))

suppressPackageStartupMessages({
  library(here)
  library(arrow)        # read_parquet
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(splines)      # ns() for vent_day spline in A3
  library(lme4)         # glmer -- A3, A4, A5, A6
  library(glmmTMB)      # truncated_nbinom2 (A5 LOS), nbinom2 (A4 Part 2)
  library(survival)     # coxph counting process (A3)
  library(tidycmprsk)   # crr -- Fine-Gray (A3)
  library(broom.mixed)  # tidy model output
  library(broom)        # tidy for coxph, glm
  library(readr)        # write_csv
  library(ggplot2)
  library(patchwork)    # multi-panel figures
  library(purrr)
  library(scales)       # percent_format() for A6 figures
  library(epiR)         # epi.ccc() -- Lin's CCC (A6)
  library(blandr)       # Bland-Altman (A6)
})

cat("=== ABT-RISE Site Analysis ===\n")
cat("Site:", site_id, "\n")
cat("Setup started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- 1.1 Output folder structure ---------------------------------------------
# out_dir is provided by ABTRISE_config.R (default: ./output_to_share).

subdirs <- c(
  "diagnostics",
  "A2_criterion/models",
  "A2_criterion/tables",
  "A2_criterion/figures",
  "A3_tte_outcomes/models",
  "A3_tte_outcomes/tables",
  "A3_tte_outcomes/figures",
  "A4_VFD_outcomes/models",
  "A4_VFD_outcomes/tables",
  "A4_VFD_outcomes/figures",
  "A5_mort_outcomes/models",
  "A5_mort_outcomes/tables",
  "A5_mort_outcomes/figures",
  "A6_benchmark_outcomes/models",
  "A6_benchmark_outcomes/tables",
  "A6_benchmark_outcomes/figures"
)

for (d in subdirs) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}
cat("Output folders created under:", out_dir, "\n\n")

# --- 1.2 Export helpers ------------------------------------------------------
# All outputs prefixed with site_id automatically.

prefix_file <- function(filename) paste0(site_id, "_", filename)

export_csv <- function(df, subfolder, filename) {
  path <- file.path(out_dir, subfolder, prefix_file(filename))
  write_csv(df, path)
  cat("  Exported:", file.path(subfolder, prefix_file(filename)), "\n")
}

export_rds <- function(obj, subfolder, filename) {
  path <- file.path(out_dir, subfolder, prefix_file(filename))
  saveRDS(obj, path)
  cat("  Exported:", file.path(subfolder, prefix_file(filename)), "\n")
}

export_png <- function(plot_obj, subfolder, filename, width = 10, height = 6) {
  path <- file.path(out_dir, subfolder, prefix_file(filename))
  ggsave(path, plot = plot_obj, width = width, height = height,
         dpi = 300, bg = "white")
  cat("  Exported:", file.path(subfolder, prefix_file(filename)), "\n")
}

# --- 1.3 Shared model helpers ------------------------------------------------

# drop_single_level(): remove factor covariates that have only one level
# at this site -- used in all analysis models to prevent contrast errors
drop_single_level <- function(vars, data) {
  vars[sapply(vars, function(v) {
    if (!v %in% names(data)) return(FALSE)
    if (!is.factor(data[[v]]) && !is.character(data[[v]])) return(TRUE)
    length(unique(na.omit(data[[v]]))) > 1
  })]
}

# --- 1.4 Shared color palette ------------------------------------------------
# Paul Tol high-contrast accessible palette
# Distinguishable for deuteranopia/protanopia (red-green colorblindness)
# and in greyscale print. Used consistently across all four scripts.
#
# Trial colors (SAT / SBT):
#   SAT = blue (#0077BB), SBT = orange (#EE7733)#   Chosen for maximum contrast under all colorblindness simulations.
#
# Four-group palette (sedation figure, eligibility groups):
#   Overall = blue, SAT-eligible = cyan, SBT-eligible = orange, Either = teal
#   All from the Tol muted palette; readable side-by-side and in greyscale.

clr_sat   <- "#0077BB"   # blue   -- SAT trial color
clr_sbt   <- "#EE7733"   # orange -- SBT trial color
clr_trial <- c(SAT = clr_sat, SBT = clr_sbt)

clr_4grp <- c(
  overall  = "#0077BB",  # blue
  sat_elig = "#33BBEE",  # cyan
  sbt_elig = "#EE7733",  # orange
  either   = "#009988"   # teal
)

# JAMA_COLORS vector (referenced by ABTRISE_02_criterion_c.R).
# Mapped to Paul Tol accessible equivalents for site distribution.
JAMA_COLORS <- c(
  "#0077BB",  # [1] dark slate  -> Tol blue

  "#EE7733",  # [2] muted orange -> Tol orange
  "#33BBEE",  # [3] JAMA blue   -> Tol cyan
  "#CC3311",  # [4] muted red   -> Tol red
  "#009988",  # [5] sage green  -> Tol teal
  "#AA3377",  # [6] muted purple -> Tol purple
  "#BBBBBB"   # [7] warm gray   -> Tol grey
)

# --- 1.5 Shared ggplot theme -------------------------------------------------
# Apply to all figures via + theme_abtrise(). Override per-plot as needed.

theme_abtrise <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position    = "bottom",
      legend.text        = element_text(size = 9),
      plot.caption       = element_text(size = 8, color = "gray40", hjust = 0),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.major.y = element_line(color = "grey92"),
      plot.title         = element_text(face = "bold", size = base_size),
      strip.text         = element_text(size = base_size, face = "bold")
    )
}

# --- 1.6 Read data files -----------------------------------------------------

cat("-- Section 1: Loading data files\n")

#PP: one row per patient, per day
df_pp_raw   <- read_parquet(file.path(data_dir, "file1_person_period.parquet"))
#hosp: one row per patient (current hospitalization)
df_hosp_raw <- read_parquet(file.path(data_dir, "file2_hospitalization_level.parquet"))

cat("File 1 (person-period):   ", nrow(df_pp_raw), "rows,",
    n_distinct(df_pp_raw$hospitalization_id), "hospitalizations\n")
cat("File 2 (hospitalization): ", nrow(df_hosp_raw), "rows\n")
cat("Hospitals in File 1:      ", n_distinct(df_pp_raw$hospital_id), "\n")
cat("Hospitals in File 2:      ", n_distinct(df_hosp_raw$hospital_id), "\n\n")

# --- 1.7 Detect single-hospital site -----------------------------------------

n_hospitals     <- n_distinct(df_pp_raw$hospital_id)
single_hospital <- (n_hospitals == 1)

if (single_hospital) {
  cat("NOTE: Single-hospital site detected.",
      "Random intercepts will be dropped; fixed-effects models only.\n\n")
}
# Random intercept term -- suppressed for single-hospital sites
# Used in all analysis model formulas
re_hosp <- if (single_hospital) "" else "(1 | hospital_id)"

# =============================================================================
# SECTION 2: DATA PREPARATION
# =============================================================================

cat("-- Section 2: Data preparation\n\n")
# --- 2.1 Population consistency check ----------------------------------------
cat("-- 2.1 Population consistency check\n")

ids_pp   <- unique(df_pp_raw$hospitalization_id)
ids_hosp <- unique(df_hosp_raw$hospitalization_id)

n_pp_not_hosp <- sum(!ids_pp %in% ids_hosp)
n_hosp_not_pp <- sum(!ids_hosp %in% ids_pp)
ids_zero_vent <- ids_hosp[!ids_hosp %in% ids_pp]

cat("Hospitalizations in File 1 not in File 2:", n_pp_not_hosp, "\n")
cat("Hospitalizations in File 2 not in File 1:", n_hosp_not_pp,
    "(zero algorithm vent days -- excluded)\n")

if (n_pp_not_hosp > 0) {
  warning("UNEXPECTED: ", n_pp_not_hosp,
          " hospitalizations in File 1 not found in File 2. ",
          "Contact coordinating center before proceeding.")
}

zero_vent_profile <- df_hosp_raw %>%
  filter(hospitalization_id %in% ids_zero_vent) %>%
  summarise(
    n                = n(),
    median_icu_los   = median(ICU_LOS,    na.rm = TRUE),
    pct_death        = round(mean(death_flag == 1, na.rm = TRUE) * 100, 1),
    median_vfd       = median(VFD_28,     na.rm = TRUE),
    pct_vfd28        = round(mean(VFD_28 == 28,   na.rm = TRUE) * 100, 1),
    median_vent_days = median(n_vent_days, na.rm = TRUE)
  )

cat("\nZero-vent-day patient profile (excluded):\n")
print(as.data.frame(zero_vent_profile))
cat("\n")

# --- 2.2 Apply zero-vent-day exclusion to both files -------------------------

cat("-- 2.2 Applying zero-vent-day exclusion\n")

df_pp_raw   <- df_pp_raw   %>% filter(hospitalization_id %in% ids_pp)
df_hosp_raw <- df_hosp_raw %>% filter(hospitalization_id %in% ids_pp)

cat("File 1 after zero-vent-day exclusion:", nrow(df_pp_raw), "rows,",
    n_distinct(df_pp_raw$hospitalization_id), "hospitalizations\n")
cat("File 2 after zero-vent-day exclusion:", nrow(df_hosp_raw),
    "hospitalizations\n\n")

export_csv(zero_vent_profile, "diagnostics", "zero_vent_day_profile.csv")

# --- 2.3 Exclusion waterfall -------------------------------------------------

waterfall <- tibble(
  step        = character(),
  file        = character(),
  n_remaining = integer(),
  n_excluded  = integer(),
  reason      = character()
)

log_step <- function(wf, step, file, n_remaining, n_excluded, reason) {
  bind_rows(wf, tibble(step, file, n_remaining, n_excluded, reason))
}

waterfall <- log_step(waterfall, "1_raw_load",
                      "File1", nrow(df_pp_raw), 0L, "Raw person-period rows loaded")
waterfall <- log_step(waterfall, "1_raw_load",
                      "File2", as.integer(nrow(df_hosp_raw) + n_hosp_not_pp),
                      as.integer(n_hosp_not_pp),
                      "Raw hospitalizations loaded; zero-vent-day patients identified")
waterfall <- log_step(waterfall, "2_zero_vent_exclusion",
                      "File2", nrow(df_hosp_raw), as.integer(n_hosp_not_pp),
                      "Excluded hospitalizations with zero algorithm vent days (fast-extubation survivors, median VFD=28)")

# --- 2.3 Impossible values check and exclusion --------------------------------
# Flags and removes rows where key variables contain TRUE data coding errors --
# values that are biologically or definitionally impossible regardless of cohort.
#
# IMPORTANT DISTINCTION -- two categories intentionally separated here:
#
#   (A) TRUE IMPOSSIBLE VALUES (handled here):
#       Values that cannot exist in any real patient. SOFA > 24, FiO2 > 1.0,
#       age < 18 in an adult ICU -- these are upstream coding errors.
#       Rows are excluded and logged in the waterfall as data errors.
#
#   (B) COHORT WINDOW / INCLUSION CRITERIA (handled elsewhere, NOT here):
#       Values that are real but outside study scope. These get their own
#       dedicated handling steps with appropriate clinical rationale:
#         - vent_day > 28: handled by filter(vent_day <= 28) in Section 2.5
#         - ICU_LOS upper cap: no cap applied -- LOS = 400 is a real patient,
#           just not flagged as an error. Lower bound (>= 1) enforced in
#           Analysis 5 dataset build as a structural cohort constraint.
#         - time_to_extubation > 28: capped at 28 in Fine-Gray construction
#           (Section 3.4 of ABTRISE_345_outcomes.R), not a data error
#         - days_to_death > 28: real patients who die after follow-up window;
#           handled as censoring in survival models, not dropped here
#         - VFD_28 > 28: definitionally impossible (kept below -- this IS
#           a coding error, not a cohort criterion)

cat("-- 2.3 Impossible values check\n")

# Only true coding errors -- variables with hard biological/definitional bounds
impossible_ranges <- list(
  # File 1 -- person-period
  SOFA_prior     = c(0, 24),    # SOFA score: 0-24 by definition
  FiO2_prior     = c(0.21, 1.0),# Fraction inspired O2: room air to 100%
  PEEP_prior     = c(0, 30),    # Clinical PEEP range; > 30 is an error
  sedation_prior = c(0, 1),     # Binary flag: only 0 or 1 valid
  # File 2 -- hospitalization
  age            = c(18, 120),  # Adult ICU cohort; < 18 = wrong patient
  CCI            = c(0, 37),    # Charlson max theoretical score
  SOFA_mean      = c(0, 24),    # Same bound as SOFA_prior
  FiO2_mean      = c(0.21, 1.0),
  PEEP_mean      = c(0, 30),
  VFD_28         = c(0, 28)     # Cannot exceed 28 by definition; > 28 = derivation error
  # NOT included (cohort window variables -- see note above):
  #   vent_day, ICU_LOS, time_to_extubation, days_to_death
)

pp_range_vars   <- intersect(names(impossible_ranges), names(df_pp_raw))
hosp_range_vars <- intersect(names(impossible_ranges), names(df_hosp_raw))

impossible_pp <- purrr::map_dfr(pp_range_vars, function(v) {
  rng    <- impossible_ranges[[v]]
  n_flag <- sum(df_pp_raw[[v]] < rng[1] | df_pp_raw[[v]] > rng[2], na.rm = TRUE)
  tibble(file = "File1", variable = v,
         min_allowed = rng[1], max_allowed = rng[2],
         obs_min = round(min(df_pp_raw[[v]], na.rm = TRUE), 3),
         obs_max = round(max(df_pp_raw[[v]], na.rm = TRUE), 3),
         n_flagged = n_flag, flag = if_else(n_flag > 0, "FLAG", "ok"))
})

impossible_hosp <- purrr::map_dfr(hosp_range_vars, function(v) {
  rng    <- impossible_ranges[[v]]
  n_flag <- sum(df_hosp_raw[[v]] < rng[1] | df_hosp_raw[[v]] > rng[2], na.rm = TRUE)
  tibble(file = "File2", variable = v,
         min_allowed = rng[1], max_allowed = rng[2],
         obs_min = round(min(df_hosp_raw[[v]], na.rm = TRUE), 3),
         obs_max = round(max(df_hosp_raw[[v]], na.rm = TRUE), 3),
         n_flagged = n_flag, flag = if_else(n_flag > 0, "FLAG", "ok"))
})

impossible_all <- bind_rows(impossible_pp, impossible_hosp)
cat("Impossible values check (true coding errors only):\n")
print(as.data.frame(impossible_all))

n_flagged_pp   <- sum(impossible_pp$n_flagged)
n_flagged_hosp <- sum(impossible_hosp$n_flagged)

if (n_flagged_pp > 0 | n_flagged_hosp > 0) {
  cat("\nWARNING:", n_flagged_pp, "rows flagged in File 1 and",
      n_flagged_hosp, "in File 2.\n")
  cat("  Flagged rows contain true data coding errors and will be excluded.\n")
  cat("  See impossible_values.csv for detail. Contact CC if n_flagged is large.\n")
} else {
  cat("\nNo impossible values detected. No rows excluded on this basis.\n")
}
cat("\n")

n_pp_before   <- nrow(df_pp_raw)
n_hosp_before <- nrow(df_hosp_raw)

for (v in pp_range_vars) {
  rng <- impossible_ranges[[v]]
  df_pp_raw <- df_pp_raw %>%
    filter(is.na(.data[[v]]) | (.data[[v]] >= rng[1] & .data[[v]] <= rng[2]))
}
for (v in hosp_range_vars) {
  rng <- impossible_ranges[[v]]
  df_hosp_raw <- df_hosp_raw %>%
    filter(is.na(.data[[v]]) | (.data[[v]] >= rng[1] & .data[[v]] <= rng[2]))
}

n_pp_excl   <- n_pp_before   - nrow(df_pp_raw)
n_hosp_excl <- n_hosp_before - nrow(df_hosp_raw)

cat("Rows excluded (data coding errors) -- File 1:", n_pp_excl,
    "| File 2:", n_hosp_excl, "\n\n")

waterfall <- log_step(waterfall, "3_data_coding_error_exclusion",
                      "File1", nrow(df_pp_raw), as.integer(n_pp_excl),
                      "Excluded rows with true data coding errors (impossible variable values: SOFA > 24, FiO2 out of range, age < 18, etc.)")
waterfall <- log_step(waterfall, "3_data_coding_error_exclusion",
                      "File2", nrow(df_hosp_raw), as.integer(n_hosp_excl),
                      "Excluded rows with true data coding errors (impossible variable values: SOFA > 24, FiO2 out of range, age < 18, VFD > 28, etc.)")

export_csv(impossible_all, "diagnostics", "impossible_values.csv")

# --- 2.4 Person-period file prep (File 1) ------------------------------------
cat("-- 2.4 Person-period file preparation\n")
df_pp <- df_pp_raw %>%
  filter(vent_day >= 1, vent_day <= 28) %>%
  left_join(
    df_hosp_raw %>%
      select(hospitalization_id, death_flag, days_to_death, age, sex, CCI),
    by = "hospitalization_id"
  ) %>%
  mutate(
    died_today             = case_when(
      death_flag == 1 & days_to_death == vent_day ~ 1L, TRUE ~ 0L),
    extubated              = as.integer(extubated),
    SAT_delivered_primary  = as.integer(SAT_delivered_primary),
    SBT_delivered_2min     = as.integer(SBT_delivered_2min),
    SAT_delivered_modified = as.integer(SAT_delivered_modified),
    SBT_delivered_5min     = as.integer(SBT_delivered_5min),
    SAT_eligible           = as.integer(SAT_eligible),
    SBT_eligible           = as.integer(SBT_eligible),
    hospital_type          = as.factor(hospital_type),
    location_type          = as.factor(location_type),
    # NEE_prior: recoded binary (NA=0 -- missing vasopressor data in ICU EHR
    # almost certainly means no vasopressor administered; avoids complete case
    # loss from NEE missingness)
    NEE_prior = if_else(is.na(NEE_prior) | NEE_prior == 0, 0L, 1L),
    # Derived flag: age < 65 (mirrors df_hosp$age_u65; used in SA_age65)
    age_u65   = as.integer(age < 65)
  )

waterfall <- log_step(waterfall, "3_vent_day_filter",
                      "File1", n_distinct(df_pp$hospitalization_id),
                      n_distinct(df_pp_raw$hospitalization_id) - n_distinct(df_pp$hospitalization_id),
                      "Restricted to vent_day 1-28")

cat("File 1 after prep:", nrow(df_pp), "person-days,",
    n_distinct(df_pp$hospitalization_id), "hospitalizations,",
    n_distinct(df_pp$hospital_id), "hospitals\n")
cat("Extubation events:", sum(df_pp$extubated,            na.rm = TRUE), "\n")
cat("Death events:     ", sum(df_pp$died_today,           na.rm = TRUE), "\n")
cat("SAT-eligible days:", sum(df_pp$SAT_eligible == 1,    na.rm = TRUE), "\n")
cat("SBT-eligible days:", sum(df_pp$SBT_eligible == 1,    na.rm = TRUE), "\n")
cat("Median vent_day:  ", median(df_pp$vent_day,          na.rm = TRUE),
    "[", quantile(df_pp$vent_day, 0.25, na.rm = TRUE), "-",
    quantile(df_pp$vent_day, 0.75, na.rm = TRUE), "]\n\n")

# --- 2.5 Hospitalization-level file prep (File 2) ----------------------------
cat("-- 2.5 Hospitalization-level file preparation\n")
df_hosp <- df_hosp_raw %>%
  mutate(
    survivor_28d    = as.integer(death_flag == 0),
    alive_28d       = as.integer(death_flag == 0),
    extubation_flag = as.integer(extubation_flag),
    death_flag      = as.integer(death_flag),
    hospital_type   = as.factor(hospital_type),
    location_type   = as.factor(location_type),
    NEE_mean        = if_else(is.na(NEE_mean) | NEE_mean == 0, 0L, 1L),
    # Derived flag: age < 65 (used in SA sensitivity analysis and Table 1)
    age_u65         = as.integer(age < 65),
    event_fg        = factor(
      case_when(
        extubation_flag == 1 ~ "extubated",
        death_flag == 1      ~ "died",
        TRUE                 ~ "censored"
      ),
      levels = c("censored", "extubated", "died")
    )
  )

# mv_count range check -- flag implausible values (count must be >= 1 by
# cohort definition; very high values may indicate data anomaly)
mv_count_check <- df_hosp %>%
  summarise(
    mv_icu_min     = min(mv_count_in_index_icu_stay,          na.rm = TRUE),
    mv_icu_max     = max(mv_count_in_index_icu_stay,          na.rm = TRUE),
    mv_icu_n_lt1   = sum(mv_count_in_index_icu_stay < 1,      na.rm = TRUE),
    mv_hosp_min    = min(mv_count_in_whole_hospitalization,    na.rm = TRUE),
    mv_hosp_max    = max(mv_count_in_whole_hospitalization,    na.rm = TRUE),
    mv_hosp_n_lt1  = sum(mv_count_in_whole_hospitalization < 1, na.rm = TRUE)
  )
cat("MV count variable range check:\n")
print(as.data.frame(mv_count_check))
if (mv_count_check$mv_icu_n_lt1 > 0 | mv_count_check$mv_hosp_n_lt1 > 0) {
  cat("WARNING: MV count < 1 detected -- contact CC.\n")
} else {
  cat("  MV count variables look clean (all >= 1).\n")
}
cat("\n")

cat("File 2 after prep:", nrow(df_hosp), "hospitalizations\n")
cat("Extubations:      ", sum(df_hosp$extubation_flag, na.rm = TRUE), "\n")
cat("Deaths:           ", sum(df_hosp$death_flag,      na.rm = TRUE), "\n")
cat("Survivors at 28d: ", sum(df_hosp$survivor_28d,    na.rm = TRUE), "\n\n")


# =============================================================================
# SEDATION OVERLAP DIAGNOSTICS
# =============================================================================
# --- 1. Hospitalization level: how many agents ever used (from df_hosp) ------
# Uses *_mean > 0 as proxy for "agent used on >= 1 vent day"

cat("-- Sedation overlap: hospitalization level\n")
sed_agents <- c("propofol", "fentanyl", "hydromorphone",
                "morphine", "lorazepam", "midazolam", "nmb")

hosp_sed <- df_hosp %>%
  mutate(
    n_agents_ever = rowSums(across(
      paste0(sed_agents, "_mean"),
      ~ . > 0
    ), na.rm = TRUE)
  )

cat("Number of distinct sedation agents used during hospitalization:\n")
print(table(hosp_sed$n_agents_ever, dnn = "n_agents"))
cat("\n")

cat("Patients with >1 agent ever used:\n")
cat(" n =", sum(hosp_sed$n_agents_ever > 1),
    "(", round(mean(hosp_sed$n_agents_ever > 1) * 100, 1), "%)\n\n")

# --- 2. Patient-day level: how many agents on same vent day (from df_pp) -----
# Uses *_prior_flag variables (binary: agent given on that day)
cat("-- Sedation overlap: patient-day level\n")

day_agents <- c("propofol_prior_flag", "fentanyl_prior_flag",
                "hydromorphone_prior_flag", "morphine_prior_flag",
                "lorazepam_prior_flag", "midazolam_prior_flag",
                "nmb_prior_flag")

pp_sed <- df_pp %>%
  mutate(
    n_agents_day = rowSums(across(all_of(day_agents), ~ . == 1), na.rm = TRUE)
  )

cat("Number of distinct sedation agents on same patient-day:\n")
print(table(pp_sed$n_agents_day, dnn = "n_agents_same_day"))
cat("\n")

cat("Patient-days with >1 agent:\n")
cat(" n =", sum(pp_sed$n_agents_day > 1),
    "(", round(mean(pp_sed$n_agents_day > 1) * 100, 1), "% of all vent days)\n\n")

# Among eligible days only
pp_sed_elig <- pp_sed %>% filter(SAT_eligible == 1 | SBT_eligible == 1)
cat("Among eligible days (SAT or SBT eligible), patient-days with >1 agent:\n")
cat(" n =", sum(pp_sed_elig$n_agents_day > 1),
    "(", round(mean(pp_sed_elig$n_agents_day > 1) * 100, 1), "% of eligible days)\n\n")

# --- 3. Co-occurrence matrix (patient-day level) ------------------------------
# Which pairs of agents most commonly co-occur on the same day?

cat("-- Agent co-occurrence matrix (patient-day level, n (%) of vent days):\n")

co_mat_n <- df_pp %>%
  select(all_of(day_agents)) %>%
  mutate(across(everything(), ~ as.integer(. == 1))) %>%
  { t(.) %*% as.matrix(.) }

co_mat_pct <- round(co_mat_n / nrow(df_pp) * 100, 1)

co_mat_combined <- matrix(
  sprintf("%d (%.1f%%)", as.integer(co_mat_n), co_mat_pct), #format to n (x%)
  nrow = nrow(co_mat_n),
  dimnames = list(gsub("_prior_flag", "", day_agents),
                  gsub("_prior_flag", "", day_agents))
)

print(co_mat_combined, quote = FALSE) # FALSE suppresses quotation marks in matrix on print 
cat("(diagonal = n days agent used; off-diagonal = n days both used together)\n\n")


cat("Zero-sedation patient profile (clinical coherence check):\n")
# Do zero-sedation patients look different from sedated patients?
zero_sed_profile <- df_hosp %>%
  mutate(any_sedation = sedation_mean > 0) %>%
  group_by(any_sedation) %>%
  summarise(
    n                = n(),
    median_vent_days = median(n_vent_days, na.rm = TRUE),
    median_sofa      = median(SOFA_mean,   na.rm = TRUE),
    pct_extubated    = round(mean(extubation_flag, na.rm = TRUE) * 100, 1),
    pct_death        = round(mean(death_flag,      na.rm = TRUE) * 100, 1),
    .groups          = "drop"
  )
print(as.data.frame(zero_sed_profile))
cat("\n")

# --- Sedation diagnostics exports (tables + figures) -----------------------

cat("-- Exporting sedation overlap diagnostics\n")

# 3a. Agent count distributions -- hospitalization and day level
sed_hosp_dist <- tibble(
  level          = "hospitalization",
  n_agents       = as.integer(names(table(hosp_sed$n_agents_ever))),
  n_observations = as.integer(table(hosp_sed$n_agents_ever)),
  pct            = round(as.numeric(table(hosp_sed$n_agents_ever)) /
                           nrow(hosp_sed) * 100, 1)
)

sed_day_dist <- tibble(
  level          = "patient_day",
  n_agents       = as.integer(names(table(pp_sed$n_agents_day))),
  n_observations = as.integer(table(pp_sed$n_agents_day)),
  pct            = round(as.numeric(table(pp_sed$n_agents_day)) /
                           nrow(pp_sed) * 100, 1)
)

sed_day_elig_dist <- tibble(
  level          = "patient_day_eligible",
  n_agents       = as.integer(names(table(pp_sed_elig$n_agents_day))),
  n_observations = as.integer(table(pp_sed_elig$n_agents_day)),
  pct            = round(as.numeric(table(pp_sed_elig$n_agents_day)) /
                           nrow(pp_sed_elig) * 100, 1)
)

export_csv(
  bind_rows(sed_hosp_dist, sed_day_dist, sed_day_elig_dist),
  "diagnostics", "sedation_agent_count_distribution.csv"
)

# 3b. Co-occurrence matrix as tidy long-format CSV
agent_labels <- gsub("_prior_flag", "", day_agents)
co_tidy <- as_tibble(co_mat_n) %>%
  mutate(agent_row = agent_labels) %>%
  pivot_longer(cols = -agent_row, names_to = "agent_col", values_to = "n_days") %>%
  mutate(
    agent_col = gsub("_prior_flag", "", agent_col),
    pct_days  = round(n_days / nrow(df_pp) * 100, 2),
    n_total_days = nrow(df_pp)
  )

export_csv(co_tidy, "diagnostics", "sedation_cooccurrence_matrix.csv")

# 3c. Zero-sedation profile
export_csv(zero_sed_profile, "diagnostics", "zero_sedation_profile.csv")

# 3d. Sedation overlap figures

# Figure: co-occurrence heatmap (lower triangle only)
co_heatmap_data <- co_tidy %>%
  # Lower triangle: agent_row >= agent_col alphabetically, excluding diagonal
  filter(agent_row != agent_col) %>%
  # Keep unique pairs (lower triangle)
  rowwise() %>%
  mutate(pair = paste(sort(c(agent_row, agent_col)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

fig_sed_heatmap <- ggplot(co_tidy,
                           aes(x = agent_col, y = agent_row,
                               fill = pct_days)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct_days)),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "gray85", high = JAMA_COLORS[3],
                      name = "% of vent days",
                      labels = function(x) paste0(round(x, 1), "%")) +
  scale_x_discrete(position = "top") +
  labs(
    title    = "Sedation Agent Co-Occurrence (Patient-Day Level)",
    subtitle = paste0("n = ", format(nrow(df_pp), big.mark = ","),
                      " vent days | diagonal = single-agent prevalence"),
    x        = NULL,
    y        = NULL
  ) +
  theme_abtrise() +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 0, size = 9),
    axis.text.y       = element_text(size = 9),
    panel.grid        = element_blank(),
    panel.border      = element_blank(),
    axis.ticks        = element_blank(),
    legend.position   = "right"
  )

export_png(fig_sed_heatmap, "diagnostics", "fig_sedation_cooccurrence.png",
           width = 7, height = 5.5)

# Figure: overlap bar chart -- stacked bars by n_agents (hosp level)
fig_sed_overlap <- ggplot(sed_hosp_dist,
                           aes(x = factor(n_agents), y = pct)) +
  geom_col(fill = JAMA_COLORS[3], alpha = 0.85, width = 0.65) +
  geom_text(aes(label = paste0(n_observations, "\n(", pct, "%)")),
            vjust = -0.3, size = 3, color = "gray25") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Number of Distinct Sedation Agents Per Hospitalization",
    subtitle = paste0("n = ", format(nrow(hosp_sed), big.mark = ","),
                      " hospitalizations"),
    x        = "Number of agents used (any vent day)",
    y        = "% of hospitalizations"
  ) +
  theme_abtrise() +
  theme(panel.grid.minor = element_blank())

export_png(fig_sed_overlap, "diagnostics", "fig_sedation_overlap_counts.png",
           width = 7, height = 5)

cat("  Sedation diagnostics exported.\n\n")

cat("-- 2.6 Covariate missingness check\n")

# ---- 2.6a Full-cohort missingness (all covariates, both files) ---------------

miss_vars_pp <- c("SOFA_prior", "FiO2_prior", "PEEP_prior",
                  "sedation_prior", "NEE_prior",
                  "SAT_delivered_primary", "SBT_delivered_2min",
                  "SAT_delivered_modified", "SBT_delivered_5min",
                  "SAT_eligible", "SBT_eligible",
                  "extubated", "vent_day")

miss_vars_hosp <- c("age", "sex", "CCI",
                    "SOFA_mean", "FiO2_mean", "PEEP_mean", "sedation_mean",
                    "NEE_mean",
                    "SAT_prop_final_primary", "SBT_prop_final_2min",
                    "SAT_prop_final_modified", "SBT_prop_final_5min",
                    "alive_28d", "VFD_28", "ICU_LOS", "death_flag",
                    "extubation_flag", "time_to_extubation",
                    "mv_count_in_index_icu_stay",
                    "mv_count_in_whole_hospitalization")

# Only check variables that exist in the data
miss_vars_pp   <- intersect(miss_vars_pp,   names(df_pp))
miss_vars_hosp <- intersect(miss_vars_hosp, names(df_hosp))

miss_pp <- df_pp %>%
  summarise(across(all_of(miss_vars_pp),
                   list(n_miss   = ~ sum(is.na(.)),
                        pct_miss = ~ round(mean(is.na(.)) * 100, 2)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(), names_to = "var_fn", values_to = "value") %>%
  separate(var_fn, into = c("variable", "fn"), sep = "__") %>%
  pivot_wider(names_from = fn, values_from = value) %>%
  mutate(
    file = "File1",
    n_total = nrow(df_pp),
    note = case_when(
      variable %in% c("SOFA_prior","FiO2_prior","PEEP_prior","sedation_prior") ~
        "Day 1 NAs expected by design",
      variable == "NEE_prior" ~
        "Binary recoded (NA=0); should show 0% missing after recode",
      TRUE ~ ""
    )
  )

miss_hosp <- df_hosp %>%
  summarise(across(all_of(miss_vars_hosp),
                   list(n_miss   = ~ sum(is.na(.)),
                        pct_miss = ~ round(mean(is.na(.)) * 100, 2)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(), names_to = "var_fn", values_to = "value") %>%
  separate(var_fn, into = c("variable", "fn"), sep = "__") %>%
  pivot_wider(names_from = fn, values_from = value) %>%
  mutate(file = "File2", n_total = nrow(df_hosp), note = "")

miss_full_cohort <- bind_rows(miss_pp, miss_hosp) %>%
  select(file, variable, n_miss, pct_miss, n_total, note) %>%
  arrange(file, desc(pct_miss))

cat("Full-cohort missingness summary:\n")
print(as.data.frame(miss_full_cohort %>% filter(pct_miss > 0)))

flag_miss <- miss_full_cohort %>%
  filter(pct_miss > 10,
         !variable %in% c("SOFA_prior","FiO2_prior","PEEP_prior",
                          "sedation_prior","NEE_prior",
                          "time_to_extubation"))

if (nrow(flag_miss) > 0) {
  cat("\nWARNING: Variables > 10% missing (excluding expected NAs):\n")
  print(as.data.frame(flag_miss))
} else {
  cat("\nAll key covariates < 10% missing. Complete case analysis proceeds.\n")
}
cat("\n")

# NEE_prior recode verification
nee_miss_check <- mean(is.na(df_pp$NEE_prior)) * 100
if (nee_miss_check > 0) {
  warning("NEE_prior still has ", round(nee_miss_check, 1),
          "% missing after binary recode -- check Section 2.4 recoding.")
} else {
  cat("NEE_prior binary recode verified: 0% missing (NA coded as 0).\n\n")
}

# ---- 2.6b Analysis-level missingness and exclusion counts --------------------
# For each analysis, compute how many observations would be excluded by
# complete-case filtering on the required covariates. This lets sites preview
# analytic sample sizes before running the downstream scripts.

cat("-- 2.6b Analysis-level complete-case preview\n\n")

# Helper: count complete and excluded for a set of covariates
cc_preview <- function(data, file_label, analysis_label, vars) {
  vars_present <- intersect(vars, names(data))
  vars_absent  <- setdiff(vars, names(data))

  n_total    <- nrow(data)
  n_complete <- sum(complete.cases(data[, vars_present, drop = FALSE]))
  n_excluded <- n_total - n_complete

  # Per-variable contribution to missingness
  per_var <- purrr::map_dfr(vars_present, function(v) {
    tibble(
      analysis  = analysis_label,
      file      = file_label,
      variable  = v,
      n_missing = sum(is.na(data[[v]])),
      pct_missing = round(mean(is.na(data[[v]])) * 100, 2)
    )
  })

  summary_row <- tibble(
    analysis       = analysis_label,
    file           = file_label,
    n_total        = n_total,
    n_complete     = n_complete,
    n_excluded     = n_excluded,
    pct_excluded   = round(n_excluded / n_total * 100, 1),
    vars_checked   = paste(vars_present, collapse = ", "),
    vars_absent    = if (length(vars_absent) > 0)
                       paste(vars_absent, collapse = ", ") else ""
  )

  list(summary = summary_row, per_var = per_var)
}

# A3 (File 1 -- person-period discrete-time):
cc_a3 <- cc_preview(
  df_pp, "File1", "A3_discrete_time",
  c("SAT_delivered_primary", "SBT_delivered_2min",
    "SOFA_prior", "FiO2_prior", "PEEP_prior", "sedation_prior",
    "age", "sex", "CCI")
)

# A4 (File 2 -- hospitalization-level two-part model):
cc_a4 <- cc_preview(
  df_hosp, "File2", "A4_VFD28_twopart",
  c("SAT_prop_final_primary", "SBT_prop_final_2min",
    "age", "sex", "CCI",
    "SOFA_mean", "FiO2_mean", "PEEP_mean", "sedation_mean",
    "alive_28d", "VFD_28")
)

# A5 (File 2 -- ICU LOS + mortality):
cc_a5 <- cc_preview(
  df_hosp, "File2", "A5_LOS_mortality",
  c("SAT_prop_final_primary", "SBT_prop_final_2min",
    "age", "sex", "CCI",
    "SOFA_mean", "FiO2_mean", "PEEP_mean", "sedation_mean",
    "ICU_LOS", "death_flag")
)

# A6 SAT (File 1 -- eligible SAT days):
cc_a6_sat <- cc_preview(
  df_pp %>% filter(SAT_eligible == 1L),
  "File1", "A6_SAT_benchmarking",
  c("SAT_delivered_primary",
    "age", "sex", "CCI",
    "SOFA_prior", "FiO2_prior", "PEEP_prior", "sedation_prior")
)

# A6 SBT (File 1 -- eligible SBT days):
cc_a6_sbt <- cc_preview(
  df_pp %>% filter(SBT_eligible == 1L),
  "File1", "A6_SBT_benchmarking",
  c("SBT_delivered_2min",
    "age", "sex", "CCI",
    "SOFA_prior", "FiO2_prior", "PEEP_prior", "sedation_prior")
)

# Combine summaries
cc_summary_all <- bind_rows(
  cc_a3$summary, cc_a4$summary, cc_a5$summary,
  cc_a6_sat$summary, cc_a6_sbt$summary
)

cc_pervar_all <- bind_rows(
  cc_a3$per_var, cc_a4$per_var, cc_a5$per_var,
  cc_a6_sat$per_var, cc_a6_sbt$per_var
) %>%
  filter(n_missing > 0)

cat("Analysis-level complete-case preview:\n")
print(as.data.frame(cc_summary_all %>%
  select(analysis, file, n_total, n_complete, n_excluded, pct_excluded)))
cat("\n")

if (nrow(cc_pervar_all) > 0) {
  cat("Variables contributing to exclusion (missing > 0):\n")
  print(as.data.frame(cc_pervar_all))
  cat("\n")
}

# Flag analyses losing > 15% to complete case
cc_high_loss <- cc_summary_all %>% filter(pct_excluded > 15)
if (nrow(cc_high_loss) > 0) {
  cat("WARNING: The following analyses lose >15% to complete-case filtering:\n")
  for (i in seq_len(nrow(cc_high_loss))) {
    cat("  ", cc_high_loss$analysis[i], ": ",
        cc_high_loss$n_excluded[i], " excluded (",
        cc_high_loss$pct_excluded[i], "%)\n")
  }
  cat("  Review per-variable missingness to identify primary driver(s).\n")
  cat("  Consider imputation or contacting CC for guidance.\n\n")
} else {
  cat("All analyses < 15% complete-case loss.\n\n")
}

# Export both tables
export_csv(miss_full_cohort, "diagnostics", "missingness_full_cohort.csv")
export_csv(cc_summary_all,   "diagnostics", "missingness_analysis_level_summary.csv")
export_csv(cc_pervar_all,    "diagnostics", "missingness_analysis_level_pervar.csv")

# --- 2.8 Delivery rate descriptives ------------------------------------------

cat("-- 2.7 Delivery rates by hospital and exposure definition\n")

delivery_rates <- df_pp %>%
  group_by(hospital_id) %>%
  summarise(
    n_vent_days            = n(),
    n_SAT_eligible         = sum(SAT_eligible,          na.rm = TRUE),
    n_SAT_primary          = sum(SAT_delivered_primary,  na.rm = TRUE),
    n_SAT_modified         = sum(SAT_delivered_modified, na.rm = TRUE),
    rate_SAT_primary       = round(n_SAT_primary  / n_SAT_eligible * 100, 1),
    rate_SAT_modified      = round(n_SAT_modified / n_SAT_eligible * 100, 1),
    n_SBT_eligible         = sum(SBT_eligible,           na.rm = TRUE),
    n_SBT_2min             = sum(SBT_delivered_2min,     na.rm = TRUE),
    n_SBT_5min             = sum(SBT_delivered_5min,     na.rm = TRUE),
    rate_SBT_2min          = round(n_SBT_2min / n_SBT_eligible * 100, 1),
    rate_SBT_5min          = round(n_SBT_5min / n_SBT_eligible * 100, 1),
    SBT_2min_to_5min_ratio = round(rate_SBT_2min / rate_SBT_5min, 2),
    .groups = "drop"
  )

cat("Delivery rates by hospital:\n")
print(as.data.frame(delivery_rates))
cat("\n")

delivery_by_day <- df_pp %>%
  group_by(vent_day) %>%
  summarise(
    n_days            = n(),
    rate_SAT_primary  = round(mean(SAT_delivered_primary  == 1 &
                                     SAT_eligible == 1, na.rm = TRUE) * 100, 1),
    rate_SAT_modified = round(mean(SAT_delivered_modified == 1 &
                                     SAT_eligible == 1, na.rm = TRUE) * 100, 1),
    rate_SBT_2min     = round(mean(SBT_delivered_2min == 1 &
                                     SBT_eligible == 1, na.rm = TRUE) * 100, 1),
    rate_SBT_5min     = round(mean(SBT_delivered_5min == 1 &
                                     SBT_eligible == 1, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

# --- 2.9 Shared covariate lists ----------------------------------------------
# Used by all analysis scripts via source()

covariates_baseline <- c("age", "sex", "CCI")

covariates_tv       <- c("SOFA_prior", "FiO2_prior", "PEEP_prior",
                         "sedation_prior", "NEE_prior")

covariates_mean     <- c("SOFA_mean", "FiO2_mean", "PEEP_mean",
                         "sedation_mean")

# A6-specific covariate list (time-varying + baseline, for vent-day model)
covariates_a6 <- c("age", "sex", "CCI",
                   "SOFA_prior", "FiO2_prior", "PEEP_prior",
                   "sedation_prior", "NEE_prior")

cat("Covariate lists:\n")
cat("  Baseline:     ", paste(covariates_baseline, collapse = ", "), "\n")
cat("  Time-varying: ", paste(covariates_tv,       collapse = ", "), "\n")
cat("  Episode mean: ", paste(covariates_mean,      collapse = ", "), "\n")
cat("  A6 (vent-day):", paste(covariates_a6,        collapse = ", "), "\n")
cat("  RE term:      ",
    if (single_hospital) "NONE (single hospital)" else re_hosp, "\n\n")

# --- 2.10 Flowsheet gate detection (Analysis 2 and A6 CCC) ------------------
# Derived from data -- sites do not configure manually.
# site_has_flowsheet_sat / _sbt = TRUE if at least one positive flowsheet
# event (== 1) exists on eligible-day rows. A column of all 0s means the
# site does not document that trial in their flowsheet -- all-zero reference
# data produces spurious accuracy/specificity with zero sensitivity, so we
# require evidence that the site actually records the trial.

cat("-- 2.8 Flowsheet gate detection\n")

site_has_flowsheet_sat <- df_pp %>%
  filter(SAT_eligible == 1) %>%
  pull(flowsheet_SAT) %>%
  { any(. == 1L, na.rm = TRUE) }

site_has_flowsheet_sbt <- df_pp %>%
  filter(SBT_eligible == 1) %>%
  pull(flowsheet_SBT) %>%
  { any(. == 1L, na.rm = TRUE) }

cat("site_has_flowsheet_sat:", site_has_flowsheet_sat, "\n")
cat("site_has_flowsheet_sbt:", site_has_flowsheet_sbt, "\n")

if (!site_has_flowsheet_sat & !site_has_flowsheet_sbt) {
  cat("  NOTE: No flowsheet data -- Analysis 2 and A6 CCC section will skip.\n")
} else {
  if (site_has_flowsheet_sat)
    cat("  SAT flowsheet present -- Analysis 2 SAT will run.\n")
  if (site_has_flowsheet_sbt)
    cat("  SBT flowsheet present -- Analysis 2 SBT will run.\n")
}
cat("\n")

# --- 2.10 Table 1: Analytic Cohort Characteristics ---------------------------
# Four-column structure: Overall | SAT-eligible | SBT-eligible | Either-eligible
# Eligibility defined at hospitalization level: >= 1 eligible vent day
# Overall = all df_hosp (all hospitalizations with >= 1 algorithm vent day)
# Either  = any_SAT_eligible | any_SBT_eligible
# Continuous variables: median (IQR)
# Categorical variables: n (%)
# Sedation variables from df_hosp: episode-level means (proportion of vent
#   days with each agent delivered)
# Outputs:
#   outputs/diagnostics/table1_cohort_characteristics.csv
#   outputs/diagnostics/table1_cohort_characteristics.html -- formatted for review

cat("-- 2.10 Table 1: Analytic cohort characteristics\n")

# --- Build hospitalization-level eligibility flags ---------------------------
# Join episode-level eligibility summary from df_pp to df_hosp
# any_SAT_eligible / any_SBT_eligible = TRUE if >= 1 eligible vent day
# any_either_eligible = union of the two (SAT OR SBT eligible)

hosp_elig <- df_pp %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_SAT_eligible    = any(SAT_eligible == 1, na.rm = TRUE),
    any_SBT_eligible    = any(SBT_eligible == 1, na.rm = TRUE),
    n_SAT_eligible_days = sum(SAT_eligible == 1, na.rm = TRUE),
    n_SBT_eligible_days = sum(SBT_eligible == 1, na.rm = TRUE),
    .groups = "drop"
  )

df_t1 <- df_hosp %>%
  left_join(hosp_elig, by = "hospitalization_id") %>%
  mutate(
    any_SAT_eligible    = replace_na(any_SAT_eligible, FALSE),
    any_SBT_eligible    = replace_na(any_SBT_eligible, FALSE),
    any_either_eligible = any_SAT_eligible | any_SBT_eligible
  )

# Column n sizes
n_overall   <- nrow(df_t1)
n_sat_elig  <- sum(df_t1$any_SAT_eligible)
n_sbt_elig  <- sum(df_t1$any_SBT_eligible)
n_either    <- sum(df_t1$any_either_eligible)
n_neither   <- sum(!df_t1$any_SAT_eligible & !df_t1$any_SBT_eligible)
n_both      <- sum(df_t1$any_SAT_eligible  & df_t1$any_SBT_eligible)
n_u65       <- sum(df_t1$age_u65 == 1, na.rm = TRUE)

cat("Table 1 cohort sizes:\n")
cat("  Overall:        ", n_overall, "\n")
cat("  SAT-eligible:   ", n_sat_elig, "\n")
cat("  SBT-eligible:   ", n_sbt_elig, "\n")
cat("  Either-eligible:", n_either, "\n")
cat("  Both-eligible:  ", n_both, "\n")
cat("  Neither:        ", n_neither, "\n")
cat("  Age < 65:       ", n_u65, "(", round(n_u65/n_overall*100,1), "% of cohort)\n\n")

# =============================================================================
# SEDATION DISTRIBUTION FIGURE -- Cleveland Dot Plot
# % of patients who received each agent >= 1 vent day
# Four eligibility groups: Overall | SAT-eligible | SBT-eligible | Either-eligible
# NMB in separate panel
# =============================================================================

library(forcats)

# --- Build agent prevalence by eligibility group -----------------------------

sed_prev <- bind_rows(
  df_hosp %>%
    summarise(across(paste0(sed_agents, "_mean"), ~ round(mean(. > 0, na.rm = TRUE) * 100, 1))) %>%
    mutate(group = paste0("Overall (N=", n_overall, ")")),

  df_t1 %>% filter(any_SAT_eligible) %>%
    summarise(across(paste0(sed_agents, "_mean"), ~ round(mean(. > 0, na.rm = TRUE) * 100, 1))) %>%
    mutate(group = paste0("SAT-Eligible (N=", n_sat_elig, ")")),

  df_t1 %>% filter(any_SBT_eligible) %>%
    summarise(across(paste0(sed_agents, "_mean"), ~ round(mean(. > 0, na.rm = TRUE) * 100, 1))) %>%
    mutate(group = paste0("SBT-Eligible (N=", n_sbt_elig, ")")),

  df_t1 %>% filter(any_either_eligible) %>%
    summarise(across(paste0(sed_agents, "_mean"), ~ round(mean(. > 0, na.rm = TRUE) * 100, 1))) %>%
    mutate(group = paste0("Either-Eligible (N=", n_either, ")"))
) %>%
  pivot_longer(
    cols      = paste0(sed_agents, "_mean"),
    names_to  = "agent",
    values_to = "pct_patients"
  ) %>%
  mutate(
    agent = gsub("_mean", "", agent),
    agent = recode(agent,
                   propofol      = "Propofol",
                   fentanyl      = "Fentanyl",
                   hydromorphone = "Hydromorphone",
                   morphine      = "Morphine",
                   lorazepam     = "Lorazepam",
                   midazolam     = "Midazolam",
                   nmb           = "NMB"
    ),
    group = factor(group, levels = c(
      paste0("Overall (N=", n_overall, ")"),
      paste0("SAT-Eligible (N=", n_sat_elig, ")"),
      paste0("SBT-Eligible (N=", n_sbt_elig, ")"),
      paste0("Either-Eligible (N=", n_either, ")")
    ))
  )

# --- Also add "No Sedation" row ----------------------------------------------

no_sed_prev <- bind_rows(
  tibble(group = paste0("Overall (N=", n_overall, ")"),
         pct_patients = round(mean(df_hosp$sedation_mean == 0, na.rm = TRUE) * 100, 1)),
  tibble(group = paste0("SAT-Eligible (N=", n_sat_elig, ")"),
         pct_patients = round(mean(df_t1$sedation_mean[df_t1$any_SAT_eligible] == 0, na.rm = TRUE) * 100, 1)),
  tibble(group = paste0("SBT-Eligible (N=", n_sbt_elig, ")"),
         pct_patients = round(mean(df_t1$sedation_mean[df_t1$any_SBT_eligible] == 0, na.rm = TRUE) * 100, 1)),
  tibble(group = paste0("Either-Eligible (N=", n_either, ")"),
         pct_patients = round(mean(df_t1$sedation_mean[df_t1$any_either_eligible] == 0, na.rm = TRUE) * 100, 1))
) %>%
  mutate(
    agent = "No Sedation",
    group = factor(group, levels = levels(sed_prev$group))
  )

# --- Split into sedation panel and NMB panel ---------------------------------

sed_plot_data <- sed_prev %>%
  filter(agent != "NMB") %>%
  bind_rows(no_sed_prev) %>%
  mutate(
    # Order agents: No Sedation at top, then by overall prevalence descending
    agent = fct_reorder(agent, pct_patients, .fun = max)
  )

nmb_plot_data <- sed_prev %>%
  filter(agent == "NMB")

# --- Color palette -----------------------------------------------------------
# Use shared clr_4grp palette defined in Section 1.4

group_colors <- c(
  clr_4grp["overall"],   # blue
  clr_4grp["sat_elig"],  # cyan
  clr_4grp["sbt_elig"],  # orange
  clr_4grp["either"]     # teal
)
names(group_colors) <- levels(sed_prev$group)

# --- Panel A: Sedation agents + No Sedation ----------------------------------

p_sed <- ggplot(sed_plot_data,
                aes(x = pct_patients, y = agent, color = group)) +
  geom_line(aes(group = agent), color = "grey80", linewidth = 0.6) +
  geom_point(size = 3.5) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title   = "A  Sedation Agent Use",
    x       = "Patients receiving agent (% of group)",
    y       = NULL
  ) +
  theme_abtrise() +
  theme(axis.text.y = element_text(size = 10)) +
  guides(color = guide_legend(nrow = 2))

# --- Panel B: NMB ------------------------------------------------------------

p_nmb <- ggplot(nmb_plot_data,
                aes(x = pct_patients, y = agent, color = group)) +
  geom_line(aes(group = agent), color = "grey80", linewidth = 0.6) +
  geom_point(size = 3.5) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "B  Neuromuscular Blockade (NMB)",
    x     = "Patients receiving NMB (% of group)",
    y     = NULL
  ) +
  theme_abtrise() +
  theme(legend.position = "none",
        axis.text.y     = element_text(size = 10))

# --- Combine with patchwork --------------------------------------------------
# Panel A taller since it has more rows; B is a slim strip

fig_sedation <- p_sed / p_nmb +
  plot_layout(heights = c(6, 1), guides = "collect") &
  theme(legend.position = "bottom")

print(fig_sedation)

export_png(fig_sedation, "diagnostics", "fig_sedation_distribution.png",
           width = 7, height = 7)

cat("Sedation distribution figure exported.\n\n")

# Export sedation prevalence data for CC aggregation
export_csv(
  bind_rows(sed_prev, no_sed_prev) %>%
    mutate(
      group_type = case_when(
        str_detect(as.character(group), "^Overall")         ~ "Overall",
        str_detect(as.character(group), "^SAT-Eligible")    ~ "SAT-Eligible",
        str_detect(as.character(group), "^SBT-Eligible")    ~ "SBT-Eligible",
        str_detect(as.character(group), "^Either-Eligible") ~ "Either-Eligible"
      ),
      n_patients = as.integer(str_extract(as.character(group), "(?<=N=)\\d+"))
    ) %>%
    select(group_type, n_patients, agent, pct_patients),
  "diagnostics", "sedation_distribution_fig.csv"
)

# --- Helper functions --------------------------------------------------------

# Median (IQR) for a continuous variable -- returns formatted string
med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  paste0(round(median(x), 1),
         " (", round(quantile(x, 0.25), 1),
         "-", round(quantile(x, 0.75), 1), ")")
}

# n (%) for a binary/logical variable -- returns formatted string
n_pct <- function(x, total) {
  n <- sum(x == 1 | x == TRUE, na.rm = TRUE)
  paste0(n, " (", round(n / total * 100, 1), "%)")
}

# n (%) for a categorical variable level
n_pct_cat <- function(x, level, total) {
  n <- sum(x == level, na.rm = TRUE)
  paste0(n, " (", round(n / total * 100, 1), "%)")
}

# Missing count
n_miss <- function(x) sum(is.na(x))

# --- Build table rows --------------------------------------------------------
# Each row: variable label | overall | sat_elig | sbt_elig

build_row <- function(label, overall_val, sat_val, sbt_val, either_val,
                      indent = FALSE, bold_label = FALSE) {
  tibble(
    variable      = if (indent) paste0("  ", label) else label,
    overall       = overall_val,
    sat_eligible  = sat_val,
    sbt_eligible  = sbt_val,
    either_eligible = either_val,
    bold          = bold_label
  )
}

d_all    <- df_t1
d_sat    <- df_t1 %>% filter(any_SAT_eligible)
d_sbt    <- df_t1 %>% filter(any_SBT_eligible)
d_either <- df_t1 %>% filter(any_either_eligible)

build_t1_rows <- function(d_all, d_sat, d_sbt, d_either,
                           n_overall, n_sat_elig, n_sbt_elig, n_either) {
  bind_rows(

  # --- HEADER ROW ---
  build_row(
    "N (hospitalizations)",
    as.character(n_overall),
    as.character(n_sat_elig),
    as.character(n_sbt_elig),
    as.character(n_either),
    bold_label = TRUE
  ),

  # --- DEMOGRAPHICS ---
  build_row("DEMOGRAPHICS", "", "", "", "", bold_label = TRUE),

  build_row("Age, median (IQR)",
            med_iqr(d_all$age), med_iqr(d_sat$age), med_iqr(d_sbt$age),
            med_iqr(d_either$age), indent = TRUE),

  build_row("Male sex, n (%)",
            n_pct(d_all$sex == "Male" | d_all$sex == "M" | d_all$sex == 1, n_overall),
            n_pct(d_sat$sex == "Male" | d_sat$sex == "M" | d_sat$sex == 1, n_sat_elig),
            n_pct(d_sbt$sex == "Male" | d_sbt$sex == "M" | d_sbt$sex == 1, n_sbt_elig),
            n_pct(d_either$sex == "Male" | d_either$sex == "M" | d_either$sex == 1, n_either),
            indent = TRUE),

  build_row("Charlson Comorbidity Index, median (IQR)",
            med_iqr(d_all$CCI), med_iqr(d_sat$CCI), med_iqr(d_sbt$CCI),
            med_iqr(d_either$CCI), indent = TRUE),

  build_row("Age < 65 years, n (%)",
            n_pct(d_all$age_u65,    n_overall),
            n_pct(d_sat$age_u65,    n_sat_elig),
            n_pct(d_sbt$age_u65,    n_sbt_elig),
            n_pct(d_either$age_u65, n_either),
            indent = TRUE),

  # --- ICU CHARACTERISTICS ---
  build_row("ICU CHARACTERISTICS", "", "", "", "", bold_label = TRUE),

  build_row("Hospital type: Academic, n (%)",
            n_pct_cat(d_all$hospital_type, "academic", n_overall),
            n_pct_cat(d_sat$hospital_type, "academic", n_sat_elig),
            n_pct_cat(d_sbt$hospital_type, "academic", n_sbt_elig),
            n_pct_cat(d_either$hospital_type, "academic", n_either),
            indent = TRUE),

  build_row("ICU type: Medical, n (%)",
            n_pct_cat(d_all$location_type, "medical_icu", n_overall),
            n_pct_cat(d_sat$location_type, "medical_icu", n_sat_elig),
            n_pct_cat(d_sbt$location_type, "medical_icu", n_sbt_elig),
            n_pct_cat(d_either$location_type, "medical_icu", n_either),
            indent = TRUE),

  build_row("ICU type: Surgical, n (%)",
            n_pct_cat(d_all$location_type, "surgical_icu", n_overall),
            n_pct_cat(d_sat$location_type, "surgical_icu", n_sat_elig),
            n_pct_cat(d_sbt$location_type, "surgical_icu", n_sbt_elig),
            n_pct_cat(d_either$location_type, "surgical_icu", n_either),
            indent = TRUE),

  build_row("ICU type: General, n (%)",
            n_pct_cat(d_all$location_type, "general_icu", n_overall),
            n_pct_cat(d_sat$location_type, "general_icu", n_sat_elig),
            n_pct_cat(d_sbt$location_type, "general_icu", n_sbt_elig),
            n_pct_cat(d_either$location_type, "general_icu", n_either),
            indent = TRUE),

  build_row("ICU type: Neuro, n (%)",
            n_pct_cat(d_all$location_type, "neuro_icu", n_overall),
            n_pct_cat(d_sat$location_type, "neuro_icu", n_sat_elig),
            n_pct_cat(d_sbt$location_type, "neuro_icu", n_sbt_elig),
            n_pct_cat(d_either$location_type, "neuro_icu", n_either),
            indent = TRUE),

  build_row("ICU type: Cardiac, n (%)",
            n_pct_cat(d_all$location_type, "cardiac_icu", n_overall),
            n_pct_cat(d_sat$location_type, "cardiac_icu", n_sat_elig),
            n_pct_cat(d_sbt$location_type, "cardiac_icu", n_sbt_elig),
            n_pct_cat(d_either$location_type, "cardiac_icu", n_either),
            indent = TRUE),

  # --- ILLNESS SEVERITY ---
  build_row("ILLNESS SEVERITY", "", "", "", "", bold_label = TRUE),

  build_row("Mean SOFA score, median (IQR)",
            med_iqr(d_all$SOFA_mean), med_iqr(d_sat$SOFA_mean), med_iqr(d_sbt$SOFA_mean),
            med_iqr(d_either$SOFA_mean), indent = TRUE),

  build_row("Mean FiO2, median (IQR)",
            med_iqr(d_all$FiO2_mean), med_iqr(d_sat$FiO2_mean), med_iqr(d_sbt$FiO2_mean),
            med_iqr(d_either$FiO2_mean), indent = TRUE),

  build_row("Mean PEEP (cmH2O), median (IQR)",
            med_iqr(d_all$PEEP_mean), med_iqr(d_sat$PEEP_mean), med_iqr(d_sbt$PEEP_mean),
            med_iqr(d_either$PEEP_mean), indent = TRUE),

  build_row("Vasopressor use (any day), n (%)",
            n_pct(d_all$NEE_mean > 0, n_overall),
            n_pct(d_sat$NEE_mean > 0, n_sat_elig),
            n_pct(d_sbt$NEE_mean > 0, n_sbt_elig),
            n_pct(d_either$NEE_mean > 0, n_either),
            indent = TRUE),

  # --- VENTILATION ---
  build_row("MECHANICAL VENTILATION", "", "", "", "", bold_label = TRUE),

  build_row("Vent days (index episode), median (IQR)",
            med_iqr(d_all$n_vent_days), med_iqr(d_sat$n_vent_days), med_iqr(d_sbt$n_vent_days),
            med_iqr(d_either$n_vent_days), indent = TRUE),

  build_row("MV episodes (index ICU stay), median (IQR)",
            med_iqr(d_all$mv_count_in_index_icu_stay),
            med_iqr(d_sat$mv_count_in_index_icu_stay),
            med_iqr(d_sbt$mv_count_in_index_icu_stay),
            med_iqr(d_either$mv_count_in_index_icu_stay),
            indent = TRUE),

  build_row("MV episodes (whole hospitalization), median (IQR)",
            med_iqr(d_all$mv_count_in_whole_hospitalization),
            med_iqr(d_sat$mv_count_in_whole_hospitalization),
            med_iqr(d_sbt$mv_count_in_whole_hospitalization),
            med_iqr(d_either$mv_count_in_whole_hospitalization),
            indent = TRUE),

  build_row("SAT-eligible days, median (IQR)",
            med_iqr(d_all$n_SAT_eligible_days),
            med_iqr(d_sat$n_SAT_eligible_days),
            med_iqr(d_sbt$n_SAT_eligible_days),
            med_iqr(d_either$n_SAT_eligible_days),
            indent = TRUE),

  build_row("SBT-eligible days, median (IQR)",
            med_iqr(d_all$n_SBT_eligible_days),
            med_iqr(d_sat$n_SBT_eligible_days),
            med_iqr(d_sbt$n_SBT_eligible_days),
            med_iqr(d_either$n_SBT_eligible_days),
            indent = TRUE),

  # --- SEDATION ---
  build_row("SEDATION", "", "", "", "", bold_label = TRUE),

  build_row("Any sedation (proportion of vent days), median (IQR)",
            med_iqr(d_all$sedation_mean), med_iqr(d_sat$sedation_mean),
            med_iqr(d_sbt$sedation_mean), med_iqr(d_either$sedation_mean),
            indent = TRUE),

  build_row("Propofol (proportion of vent days), median (IQR)",
            med_iqr(d_all$propofol_mean), med_iqr(d_sat$propofol_mean),
            med_iqr(d_sbt$propofol_mean), med_iqr(d_either$propofol_mean),
            indent = TRUE),

  build_row("Midazolam (proportion of vent days), median (IQR)",
            med_iqr(d_all$midazolam_mean), med_iqr(d_sat$midazolam_mean),
            med_iqr(d_sbt$midazolam_mean), med_iqr(d_either$midazolam_mean),
            indent = TRUE),

  build_row("Lorazepam (proportion of vent days), median (IQR)",
            med_iqr(d_all$lorazepam_mean), med_iqr(d_sat$lorazepam_mean),
            med_iqr(d_sbt$lorazepam_mean), med_iqr(d_either$lorazepam_mean),
            indent = TRUE),

  build_row("Fentanyl (proportion of vent days), median (IQR)",
            med_iqr(d_all$fentanyl_mean), med_iqr(d_sat$fentanyl_mean),
            med_iqr(d_sbt$fentanyl_mean), med_iqr(d_either$fentanyl_mean),
            indent = TRUE),

  build_row("Hydromorphone (proportion of vent days), median (IQR)",
            med_iqr(d_all$hydromorphone_mean), med_iqr(d_sat$hydromorphone_mean),
            med_iqr(d_sbt$hydromorphone_mean), med_iqr(d_either$hydromorphone_mean),
            indent = TRUE),

  build_row("Morphine (proportion of vent days), median (IQR)",
            med_iqr(d_all$morphine_mean), med_iqr(d_sat$morphine_mean),
            med_iqr(d_sbt$morphine_mean), med_iqr(d_either$morphine_mean),
            indent = TRUE),

  build_row("Neuromuscular blockade (proportion of vent days), median (IQR)",
            med_iqr(d_all$nmb_mean), med_iqr(d_sat$nmb_mean),
            med_iqr(d_sbt$nmb_mean), med_iqr(d_either$nmb_mean),
            indent = TRUE),

  # --- SAT/SBT DELIVERY ---
  build_row("SAT/SBT DELIVERY", "", "", "", "", bold_label = TRUE),

  build_row("SAT delivery (primary, proportion of eligible days), median (IQR)",
            med_iqr(d_all$SAT_prop_final_primary),
            med_iqr(d_sat$SAT_prop_final_primary),
            med_iqr(d_sbt$SAT_prop_final_primary),
            med_iqr(d_either$SAT_prop_final_primary),
            indent = TRUE),

  build_row("SBT delivery (2-min primary, proportion of eligible days), median (IQR)",
            med_iqr(d_all$SBT_prop_final_2min),
            med_iqr(d_sat$SBT_prop_final_2min),
            med_iqr(d_sbt$SBT_prop_final_2min),
            med_iqr(d_either$SBT_prop_final_2min),
            indent = TRUE),

  # --- OUTCOMES ---
  build_row("OUTCOMES", "", "", "", "", bold_label = TRUE),

  build_row("ICU LOS (days), median (IQR)",
            med_iqr(d_all$ICU_LOS), med_iqr(d_sat$ICU_LOS), med_iqr(d_sbt$ICU_LOS),
            med_iqr(d_either$ICU_LOS), indent = TRUE),

  build_row("Ventilator-free days at 28d, median (IQR)",
            med_iqr(d_all$VFD_28), med_iqr(d_sat$VFD_28), med_iqr(d_sbt$VFD_28),
            med_iqr(d_either$VFD_28), indent = TRUE),

  build_row("Extubated within 28d, n (%)",
            n_pct(d_all$extubation_flag, n_overall),
            n_pct(d_sat$extubation_flag, n_sat_elig),
            n_pct(d_sbt$extubation_flag, n_sbt_elig),
            n_pct(d_either$extubation_flag, n_either),
            indent = TRUE),

  build_row("In-hospital death, n (%)",
            n_pct(d_all$death_flag, n_overall),
            n_pct(d_sat$death_flag, n_sat_elig),
            n_pct(d_sbt$death_flag, n_sbt_elig),
            n_pct(d_either$death_flag, n_either),
            indent = TRUE)

) %>%
    select(-bold)
}

t1_rows <- build_t1_rows(d_all, d_sat, d_sbt, d_either,
                          n_overall, n_sat_elig, n_sbt_elig, n_either)

# --- Rename columns with n in header -----------------------------------------
names(t1_rows) <- c(
  "Variable",
  paste0("Overall (N=", n_overall, ")"),
  paste0("SAT-Eligible (N=", n_sat_elig, ")"),
  paste0("SBT-Eligible (N=", n_sbt_elig, ")"),
  paste0("Either-Eligible (N=", n_either, ")")
)

cat("Table 1 built:", nrow(t1_rows), "rows\n")
print(as.data.frame(t1_rows), right = FALSE)
cat("\n")

# --- Add footnote row --------------------------------------------------------
# Append footnote as final rows for transparency
t1_footnotes <- tibble(
  Variable = c(
    "",
    paste0("* SAT-eligible: >= 1 vent day meeting SAT eligibility criteria (n=",
           n_sat_elig, "; ", round(n_sat_elig/n_overall*100,1), "% of cohort)"),
    paste0("* SBT-eligible: >= 1 vent day meeting SBT eligibility criteria (n=",
           n_sbt_elig, "; ", round(n_sbt_elig/n_overall*100,1), "% of cohort)"),
    paste0("* Either-eligible: >= 1 SAT- or SBT-eligible vent day (n=",
           n_either, "; ", round(n_either/n_overall*100,1), "% of cohort)"),
    paste0("* Both eligible: ", n_both,
           " (", round(n_both/n_overall*100,1), "%)"),
    paste0("* Neither eligible: ", n_neither,
           " (", round(n_neither/n_overall*100,1),
           "%) -- no SAT or SBT eligible vent days in index episode"),
    paste0("* Age < 65: ", n_u65,
           " (", round(n_u65/n_overall*100,1),
           "%) -- derived flag used in sensitivity analysis SA_age65"),
    "* Sedation variables: proportion of algorithm-captured vent days (days 1-28) with agent delivered",
    "* Continuous variables: median (IQR). Categorical: n (%)",
    paste0("* Site: ", site_id, " | Single hospital: ", single_hospital,
           " | Run date: ", format(Sys.Date(), "%Y-%m-%d"))
  )
)

# Pad footnote tibble to match column count
for (col in names(t1_rows)[-1]) {
  t1_footnotes[[col]] <- ""
}
t1_footnotes <- t1_footnotes[, names(t1_rows)]

t1_export <- bind_rows(t1_rows, t1_footnotes)

# --- Export CSV --------------------------------------------------------------
export_csv(t1_export, "diagnostics", "table1_cohort_characteristics.csv")

# --- SA: Table 1 split by age group (<65 vs >=65) ----------------------------
cat("-- 2.10b Table 1 by age group (sensitivity analysis)\n")

df_t1_u65  <- df_t1 %>% filter(age <  65)
df_t1_ge65 <- df_t1 %>% filter(age >= 65)

n_overall_u65   <- nrow(df_t1_u65)
n_sat_elig_u65  <- sum(df_t1_u65$any_SAT_eligible)
n_sbt_elig_u65  <- sum(df_t1_u65$any_SBT_eligible)
n_either_u65    <- sum(df_t1_u65$any_either_eligible)

n_overall_ge65  <- nrow(df_t1_ge65)
n_sat_elig_ge65 <- sum(df_t1_ge65$any_SAT_eligible)
n_sbt_elig_ge65 <- sum(df_t1_ge65$any_SBT_eligible)
n_either_ge65   <- sum(df_t1_ge65$any_either_eligible)

t1_u65 <- build_t1_rows(
  df_t1_u65, df_t1_u65 %>% filter(any_SAT_eligible),
  df_t1_u65 %>% filter(any_SBT_eligible), df_t1_u65 %>% filter(any_either_eligible),
  n_overall_u65, n_sat_elig_u65, n_sbt_elig_u65, n_either_u65
)
names(t1_u65) <- c(
  "Variable",
  paste0("u65_Overall (N=", n_overall_u65, ")"),
  paste0("u65_SAT-Eligible (N=", n_sat_elig_u65, ")"),
  paste0("u65_SBT-Eligible (N=", n_sbt_elig_u65, ")"),
  paste0("u65_Either-Eligible (N=", n_either_u65, ")")
)

t1_ge65 <- build_t1_rows(
  df_t1_ge65, df_t1_ge65 %>% filter(any_SAT_eligible),
  df_t1_ge65 %>% filter(any_SBT_eligible), df_t1_ge65 %>% filter(any_either_eligible),
  n_overall_ge65, n_sat_elig_ge65, n_sbt_elig_ge65, n_either_ge65
)
names(t1_ge65) <- c(
  "Variable",
  paste0("ge65_Overall (N=", n_overall_ge65, ")"),
  paste0("ge65_SAT-Eligible (N=", n_sat_elig_ge65, ")"),
  paste0("ge65_SBT-Eligible (N=", n_sbt_elig_ge65, ")"),
  paste0("ge65_Either-Eligible (N=", n_either_ge65, ")")
)

t1_age_split <- bind_cols(t1_u65, t1_ge65 %>% select(-Variable))
export_csv(t1_age_split, "diagnostics", "SA_age_split_table1_cohort_characteristics.csv")
cat("  Age-split Table 1 exported: <65 n=", n_overall_u65,
    "| >=65 n=", n_overall_ge65, "\n\n")

# --- Formatted flextable (HTML preview) --------------------------------------
# Requires flextable package -- skip gracefully if not available
tryCatch({
  library(flextable)
  
  # Bold helper: rows that are section headers (no indent, no data)
  header_rows <- which(
    t1_rows[[2]] == "" &
      !grepl("^  ", t1_rows[[1]]) &
      t1_rows[[1]] != paste0("N (hospitalizations)")
  )
  
  ft <- flextable(t1_rows) %>%
    bold(i = c(1, header_rows), bold = TRUE) %>%
    bold(part = "header", bold = TRUE) %>%
    bg(i = header_rows, bg = "#D6EAF8") %>%
    bg(i = 1, bg = "#1A5276") %>%
    color(i = 1, color = "white") %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    padding(padding = 3, part = "all") %>%
    hline(i = header_rows, border = officer::fp_border(color = "#AED6F1", width = 1)) %>%
    add_footer_lines(paste0(
      "SAT-eligible: >= 1 SAT-eligible vent day. ",
      "SBT-eligible: >= 1 SBT-eligible vent day. ",
      "Either-eligible: >= 1 SAT- or SBT-eligible vent day (union). ",
      "Sedation: proportion of vent days with agent delivered. ",
      "Continuous: median (IQR). Categorical: n (%). ",
      "Site: ", site_id, " | Run: ", format(Sys.Date(), "%Y-%m-%d"), "."
    )) %>%
    fontsize(i = 1, size = 7, part = "footer")
  
  # Save as HTML for browser preview
  html_path <- file.path(out_dir, "diagnostics",
                         paste0(site_id, "_table1_cohort_characteristics.html"))
  save_as_html(ft, path = html_path)
  cat("  Table 1 HTML preview saved:", html_path, "\n")
  
}, error = function(e) {
  cat("  NOTE: flextable HTML export skipped:", conditionMessage(e), "\n")
  cat("  CSV export completed successfully -- HTML requires flextable package.\n")
})

cat("Table 1 complete.\n\n")

# --- 2.11 Diagnostic exports -------------------------------------------------

cat("-- 2.9 Exporting setup diagnostics\n")

export_csv(miss_full_cohort, "diagnostics", "missingness_analytic.csv")
export_csv(waterfall,      "diagnostics", "exclusion_waterfall.csv")
export_csv(delivery_rates, "diagnostics", "delivery_rates_by_hospital.csv")
export_csv(delivery_by_day,"diagnostics", "delivery_rates_fig.csv")

# SA: Delivery rates by vent day split by age group (<65 vs >=65)
summarise_delivery_by_day <- function(df) {
  df %>%
    group_by(vent_day) %>%
    summarise(
      n_days            = n(),
      rate_SAT_primary  = round(mean(SAT_delivered_primary  == 1 &
                                       SAT_eligible == 1, na.rm = TRUE) * 100, 1),
      rate_SAT_modified = round(mean(SAT_delivered_modified == 1 &
                                       SAT_eligible == 1, na.rm = TRUE) * 100, 1),
      rate_SBT_2min     = round(mean(SBT_delivered_2min == 1 &
                                       SBT_eligible == 1, na.rm = TRUE) * 100, 1),
      rate_SBT_5min     = round(mean(SBT_delivered_5min == 1 &
                                       SBT_eligible == 1, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    )
}

delivery_by_day_age_split <- bind_rows(
  df_pp %>% filter(age <  65) %>% summarise_delivery_by_day() %>% mutate(age_group = "age_u65"),
  df_pp %>% filter(age >= 65) %>% summarise_delivery_by_day() %>% mutate(age_group = "age_ge65")
)
export_csv(delivery_by_day_age_split, "diagnostics", "SA_age_split_delivery_rates_fig.csv")

# Raw SAT/SBT delivery rates on eligible days by age group (<65 vs >=65)
# Unadjusted -- for descriptive comparison only
delivery_by_age_group <- bind_rows(
  df_pp %>%
    filter(SAT_eligible == 1) %>%
    mutate(age_group = if_else(age >= 65, "age_ge65", "age_u65")) %>%
    group_by(age_group) %>%
    summarise(
      trial           = "SAT",
      n_eligible_days = n(),
      n_delivered     = sum(SAT_delivered_primary == 1, na.rm = TRUE),
      rate_pct        = round(mean(SAT_delivered_primary == 1, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    ),
  df_pp %>%
    filter(SBT_eligible == 1) %>%
    mutate(age_group = if_else(age >= 65, "age_ge65", "age_u65")) %>%
    group_by(age_group) %>%
    summarise(
      trial           = "SBT",
      n_eligible_days = n(),
      n_delivered     = sum(SBT_delivered_2min == 1, na.rm = TRUE),
      rate_pct        = round(mean(SBT_delivered_2min == 1, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    )
) %>%
  select(trial, age_group, n_eligible_days, n_delivered, rate_pct)
export_csv(delivery_by_age_group, "diagnostics", "SA_delivery_rates_by_age_group.csv")

fig_delivery <- delivery_by_day %>%
  pivot_longer(cols = starts_with("rate_"),
               names_to = "definition", values_to = "rate") %>%
  mutate(
    trial     = if_else(str_detect(definition, "SAT"), "SAT", "SBT"),
    def_label = case_when(
      definition == "rate_SAT_primary"  ~ "SAT primary",
      definition == "rate_SAT_modified" ~ "SAT modified",
      definition == "rate_SBT_2min"     ~ "SBT 2min (primary)",
      definition == "rate_SBT_5min"     ~ "SBT 5min (sensitivity)"
    ),
    linetype  = if_else(str_detect(definition, "modified|2min"), "dashed", "solid")
  ) %>%
  ggplot(aes(x = vent_day, y = rate, color = trial, linetype = linetype)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = clr_trial) +
  scale_linetype_identity() +
  labs(
    title    = "SAT and SBT Daily Delivery Rates by Ventilator Day",
    subtitle = "Solid = primary definition | Dashed = sensitivity definition",
    x        = "Ventilator Day",
    y        = "Delivery Rate (% eligible days)",
    color    = "Trial Type",
    caption  = "Note: SBT 2min = primary definition; SBT 5min = sensitivity. Gap quantifies threshold impact."
  ) +
  theme_abtrise()

export_png(fig_delivery, "diagnostics", "fig_delivery_rates.png")

cat("\nSetup complete. Objects available for analysis scripts:\n")
cat("  df_pp, df_hosp, df_pp_raw, df_hosp_raw\n")
cat("  waterfall, log_step, n_hospitals, single_hospital, re_hosp\n")
cat("  covariates_baseline, covariates_tv, covariates_mean, covariates_a6\n")
cat("  site_has_flowsheet_sat, site_has_flowsheet_sbt\n")
cat("  export_csv, export_rds, export_png, prefix_file, drop_single_level\n")
cat("  out_dir, site_id\n")
cat("  Derived flags: df_hosp$age_u65 (age < 65), df_pp$NEE_prior (binary recode)\n")
cat("  Shared palette: clr_sat, clr_sbt, clr_trial, clr_4grp\n")
cat("  Shared theme:   theme_abtrise()\n")
cat("\nRun finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
