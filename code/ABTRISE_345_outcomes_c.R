# =============================================================================
# ABT-RISE: Site-Level Analysis Script 3 of 4
# ANALYSES 3, 4, 5 -- Construct Validity Outcomes
#
# SCRIPTS IN THIS SERIES:
#   ABTRISE_01_setup.R          -- run directly to review data quality
#   ABTRISE_02_criterion.R
#   ABTRISE_345_outcomes.R      <- YOU ARE HERE
#   ABTRISE_06_benchmarking.R
#
# HOW TO RUN:
#   Open this file and click Source (or run source("ABTRISE_345_outcomes.R"))
#   Setup runs automatically -- do NOT run ABTRISE_01_setup.R separately first.
#
# WHAT THIS SCRIPT PRODUCES:
#   outputs/A3_tte_outcomes/models/   A3_dt_primary_coefs.csv, A3_dt_primary_re_variance.csv,
#                                    A3_cox_secondary_coefs.csv, A3_fg_secondary_coefs.csv,
#                                    A3_fg_cumulative_incidence.csv, A3_sensitivity_coefs.csv,
#                                    SA_age65_A3_dt_coefs.csv, A3_fit_*.rds
#   outputs/A3_tte_outcomes/figures/  fig_A3_dt_forest.png, A3_fig_cif_curves.csv
#   outputs/A4_VFD_outcomes/models/   A4_part1_alive28d_coefs.csv, A4_part2_vfd_survivors_coefs.csv,
#                                    A4_vfd28_descriptive.csv, A4_sensitivity_coefs.csv,
#                                    SA_age65_A4_coefs.csv, A4_fit_*.rds
#   outputs/A4_VFD_outcomes/tables/   vfd28_descriptive.csv
#   outputs/A4_VFD_outcomes/figures/  fig_A4_twopart.png, A4_fig_vfd_distribution.csv
#   outputs/A5_mort_outcomes/models/  A5_icu_los_coefs.csv, A5_icu_los_overdispersion.csv,
#                                    A5_mortality_coefs.csv, A5_sensitivity_coefs.csv,
#                                    SA_age65_A5_coefs.csv, A5_fit_*.rds
#   outputs/A5_mort_outcomes/figures/ fig_A5_los_mortality.png, A5_fig_los_distribution.csv
#   outputs/diagnostics/              fig_adj_comparison.png, fig_adj_comparison_data.csv,
#                                    session_info_a345.txt
#
# ANALYSIS 3: Time to extubation (primary construct validity)
#   3.1  Discrete-time logistic regression (PRIMARY)
#   3.2  Cause-specific Cox -- secondary
#   3.3  Fine-Gray subdistribution hazard -- secondary
#   3.S  Sensitivity: alternative exposure definitions
#
# ANALYSIS 4: VFD-28 two-part model (replaces ZINB per Hajage NEJM Evidence 2025)
#   4.1  Part 1: Mixed-effects logistic (alive at 28 days)
#   4.2  Part 2: Mixed-effects NB (VFDs among survivors)
#   4.3  Descriptive VFD-28 by delivery group
#   4.S  Sensitivity: alternative exposure definitions
#
# ANALYSIS 5: ICU LOS (zero-truncated NB) and in-hospital mortality (mixed logistic)
#   5.1  ICU LOS: ZTNB
#   5.2  Mortality: mixed-effects logistic
#   5.S  Sensitivity: alternative exposure definitions
#
# SENSITIVITY ANALYSIS 6 (SA_age65): Age < 65 subgroup
#   Repeats primary A3 discrete-time, A4 two-part, and A5 LOS + mortality
#   models restricted to patients aged < 65. Rationale: construct validity
#   consistency check across age strata. Age retained as covariate within
#   the restricted cohort. Outputs: SA_age65_A3/A4/A5_coefs.csv.
#
# COORDINATING CENTER: Rush
# =============================================================================

# --- Load setup (runs Sections 0-2 automatically) ----------------------------
source(here::here("code", "ABTRISE_01_setup_c.R"))


cat("============================================================\n")
cat("ANALYSIS 3: Time to Extubation\n")
cat("============================================================\n\n")

# --- 3.1 Discrete-time dataset construction -----------------------------------
# One row per patient per vent_day up to day 28 or first event
# Extubation day included as final row (extubated = 1 allowed on same day as SAT/SBT)
# Death treated as censoring in primary discrete-time model
# (competing risk addressed directly in Fine-Gray secondary)

cat("-- 3.1 Constructing discrete-time dataset\n")

last_day <- suppressWarnings(
  df_pp %>%
    group_by(hospitalization_id) %>%
    summarise(
      first_extub = min(vent_day[extubated == 1],  na.rm = TRUE),
      first_death = min(vent_day[died_today == 1], na.rm = TRUE),
      .groups = "drop"
    )
) %>%
  mutate(
    first_extub = if_else(is.infinite(first_extub), NA_real_, first_extub),
    first_death = if_else(is.infinite(first_death), NA_real_, first_death),
    last_day    = pmin(first_extub, first_death, 28, na.rm = TRUE)
  )

df_dt <- df_pp %>%
  left_join(last_day %>% select(hospitalization_id, last_day),
            by = "hospitalization_id") %>%
  filter(vent_day <= last_day) %>%
  filter(
    !is.na(SAT_delivered_primary), !is.na(SBT_delivered_2min),
    !is.na(SOFA_prior), !is.na(FiO2_prior), !is.na(PEEP_prior),
    !is.na(sedation_prior)
    # hospital_type/location_type no longer covariates here (Thread 4
    # amendment) -- complete-case filter on them removed accordingly
  )

# Log complete case step to waterfall
waterfall <- log_step(waterfall, "5_complete_case_A3",
                      "File1", n_distinct(df_dt$hospitalization_id),
                      n_distinct(df_pp$hospitalization_id) - n_distinct(df_dt$hospitalization_id),
                      "Complete case filter: Analysis 3 discrete-time model covariates")

cat("Discrete-time dataset:    ", nrow(df_dt), "person-day rows\n")
cat("Unique hospitalizations:  ",
    n_distinct(df_dt$hospitalization_id), "\n")
cat("Extubation events:        ",
    sum(df_dt$extubated, na.rm = TRUE), "\n\n")

# --- 3.2 PRIMARY: Discrete-time logistic regression --------------------------
# glmer: extubated ~ SAT + SBT + ns(vent_day,3) + time-varying covariates
#                  + baseline covariates (age, sex, CCI; hospital_type/
#                    location_type removed -- collinear with (1|hospital_id),
#                    Thread 4 amendment) + (1|hospital_id)
# ns(vent_day, df=3) approximates flexible baseline hazard -- confirmed Q3

# Death treated as censoring (row truncation above)

cat("-- 3.2 Primary: Discrete-time logistic regression\n")

# Check for single-level factors -- drop from fixed effects if only one level
# present at this site (contrast cannot be estimated).
# drop_single_level() is defined in ABTRISE_01_setup.R and available via source().
covariates_baseline_dt <- drop_single_level(covariates_baseline, df_dt)
covariates_tv_dt       <- drop_single_level(covariates_tv,       df_dt)

if (length(covariates_baseline_dt) < length(covariates_baseline)) {
  dropped <- setdiff(covariates_baseline, covariates_baseline_dt)
  cat("NOTE: Dropping single-level covariates from A3 discrete-time model:\n")
  cat("  ", paste(dropped, collapse = ", "), "\n")
  cat("  These variables have only one level at this site and cannot\n")
  cat("  be estimated as fixed effects. Logged in site_metadata.\n\n")
}

f_dt_primary <- reformulate(
  termlabels = c(
    "SAT_delivered_primary",
    "SBT_delivered_2min",
    "ns(vent_day, df = 3)",
    covariates_tv_dt,
    covariates_baseline_dt,
    if (!single_hospital) re_hosp
  ),
  response = "extubated"
)
cat("Formula:\n"); print(f_dt_primary); cat("\n")

fit_dt_primary_warnings <- character(0)

fit_dt_primary <- withCallingHandlers(
  if (single_hospital) {
    glm(f_dt_primary, data = df_dt, family = binomial(link = "logit"))
  } else {
    glmer(f_dt_primary, data = df_dt, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 2e5)))
  },
  warning = function(w) {
    fit_dt_primary_warnings <<- c(fit_dt_primary_warnings,
                                  conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_dt_primary_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_dt_primary_warnings, collapse = "\n"), "\n")
}

tidy_dt_primary <- if (single_hospital) {
  broom::tidy(fit_dt_primary, conf.int = TRUE, exponentiate = TRUE)
} else {
  broom.mixed::tidy(fit_dt_primary, conf.int = TRUE, exponentiate = TRUE,
                    effects = "fixed")
}

cat("\nDiscrete-time model -- primary exposure terms:\n")
tidy_dt_primary %>%
  filter(str_detect(term, "SAT_delivered|SBT_delivered")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# Random effect variance
re_var_dt <- if (!single_hospital) {
  vc <- as.data.frame(VarCorr(fit_dt_primary))
  tibble(
    model        = "A3_discrete_time_primary",
    re_term      = "hospital_id",
    variance     = vc$vcov[1],
    sd           = vc$sdcor[1],
    icc          = vc$vcov[1] / (vc$vcov[1] + (pi^2 / 3)),
    n_hospitals  = n_hospitals,
    single_hosp  = FALSE,
    model_type   = "glmer",
    warnings     = paste(fit_dt_primary_warnings, collapse = "; ")
  )
} else {
  tibble(
    model        = "A3_discrete_time_primary",
    re_term      = "hospital_id",
    variance     = NA_real_,
    sd           = NA_real_,
    icc          = NA_real_,
    n_hospitals  = n_hospitals,
    single_hosp  = TRUE,
    model_type   = "glm_fallback",
    warnings     = paste(fit_dt_primary_warnings, collapse = "; ")
  )
}

cat("Random effect variance (hospital):",
    if (!single_hospital) round(re_var_dt$variance, 4) else "NA (single hospital)",
    "\n\n")

# --- 3.3 SECONDARY: Cause-specific Cox (counting process) --------------------
# Counting process: (tstart, tstop, event_cs) per eligible vent-day
# Death censored (cause-specific for extubation)
# Shared gamma frailty for hospital clustering


cat("-- 3.3 Secondary: Cause-specific Cox (counting process)\n")

df_cox <- df_dt %>%
  mutate(
    tstart   = vent_day - 1,
    tstop    = vent_day,
    event_cs = as.integer(extubated == 1 & died_today == 0)
  )

# Apply single-level factor check to Cox covariates
covariates_tv_cox       <- drop_single_level(covariates_tv,       df_cox)
covariates_baseline_cox <- drop_single_level(covariates_baseline, df_cox)

if (length(covariates_baseline_cox) < length(covariates_baseline)) {
  dropped_cox <- setdiff(covariates_baseline, covariates_baseline_cox)
  cat("NOTE: Dropping single-level covariates from A3 Cox model:",
      paste(dropped_cox, collapse = ", "), "\n")
}

f_cox <- reformulate(
  termlabels = c(
    "SAT_delivered_primary",
    "SBT_delivered_2min",
    covariates_tv_cox,
    covariates_baseline_cox,
    if (!single_hospital) "frailty(hospital_id, distribution = 'gamma')"
  ),
  response = "Surv(tstart, tstop, event_cs)"
)

cat("Formula:\n"); print(f_cox); cat("\n")

fit_cox_warnings <- character(0)

fit_cox <- withCallingHandlers(
  coxph(f_cox, data = df_cox, ties = "efron"),
  warning = function(w) {
    fit_cox_warnings <<- c(fit_cox_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_cox_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_cox_warnings, collapse = "\n"), "\n")
}

tidy_cox <- broom::tidy(fit_cox, conf.int = TRUE, exponentiate = TRUE)

cat("\nCause-specific Cox -- primary exposure terms:\n")
tidy_cox %>%
  filter(str_detect(term, "SAT_delivered|SBT_delivered")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- 3.4 SECONDARY: Fine-Gray subdistribution hazard (File 2) ----------------
# One row per patient; competing risks: extubated (1) vs died (2)
# Exposure: SAT_prop_final_primary, SBT_prop_final_2min (episode-level summary)
# Cannot take true time-varying covariates -- episode means used
# cluster() term dropped for single-hospital sites
# hospital_type/location_type removed as covariates -- collinear with
# cluster(hospital_id) term (Thread 4 amendment)
# NEE_mean excluded -- redundant with SOFA_mean

cat("-- 3.4 Secondary: Fine-Gray subdistribution hazard\n")
df_fg <- df_hosp %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean)
    # NOTE: do NOT filter on !is.na(time_to_extubation) --
    # deaths and censored patients have NA here and must be retained
    # for the competing risk structure to work correctly
  ) %>%
  mutate(
    # Build time_fg: use time_to_extubation for extubated patients,
    # days_to_death for competing event patients,
    # n_vent_days capped at 28 for censored patients
    time_fg = case_when(
      extubation_flag == 1              ~ pmin(time_to_extubation, 28),
      death_flag == 1                   ~ pmin(days_to_death, 28),
      TRUE                              ~ pmin(n_vent_days, 28)
    )
    # event_fg already constructed correctly in Section 2.5 -- do not rebuild
  )
cat("Fine-Gray dataset:        ", nrow(df_fg), "patients\n")
cat("Event distribution:\n")
print(table(df_fg$event_fg, useNA = "always"))
cat("NOTE: censored + extubated + died should sum to n patients.",
    "If died = 0, check time_to_extubation filter -- deaths have NA",
    "for time_to_extubation and must not be excluded.\n\n")

# Apply single-level factor check to Fine-Gray covariates
covars_episode_fg <- drop_single_level(
  c("age", "sex", "CCI",
    "SOFA_mean", "FiO2_mean", "PEEP_mean", "sedation_mean"),
  df_fg
)

if (length(covars_episode_fg) < 7) {
  dropped_fg <- setdiff(
    c("age","sex","CCI",
      "SOFA_mean","FiO2_mean","PEEP_mean","sedation_mean"),
    covars_episode_fg
  )
  cat("NOTE: Dropping single-level covariates from Fine-Gray model:",
      paste(dropped_fg, collapse = ", "), "\n")
}

# Build formula -- drop cluster() for single-hospital sites
fg_rhs <- c(
  "SAT_prop_final_primary",
  "SBT_prop_final_2min",
  covars_episode_fg,
  if (!single_hospital) "cluster(hospital_id)"
)

fg_formula <- as.formula(
  paste("Surv(time_fg, event_fg) ~", paste(fg_rhs, collapse = " + "))
)

cat("Formula:\n"); print(fg_formula); cat("\n")

fit_fg_warnings <- character(0)

fit_fg <- withCallingHandlers(
  crr(fg_formula, data = df_fg, failcode = "extubated"),
  warning = function(w) {
    fit_fg_warnings <<- c(fit_fg_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_fg_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_fg_warnings, collapse = "\n"), "\n")
}

tidy_fg <- broom::tidy(fit_fg, conf.int = TRUE, exponentiate = TRUE)

cat("Fine-Gray -- primary exposure terms:\n")
tidy_fg %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# Cumulative incidence curve data

cif_data <- tryCatch({
  # Extract predicted cumulative incidence at observed times
  # ntime = NULL returns all unique event times in the dataset
  tidy_cif <- tidycmprsk::tidy(fit_fg)
  as.data.frame(tidy_cif)
}, error = function(e) {
  cat("NOTE: CIF extraction via tidy() failed:", conditionMessage(e), "\n")
  cat("  Trying alternative extraction via fit_fg object slots...\n")
  # Fallback: extract from the underlying cmprsk object if accessible
  if (!is.null(fit_fg$cmprsk)) {
    as.data.frame(fit_fg$cmprsk$uftime)
  } else {
    cat("  CIF data not extractable -- exporting empty placeholder.\n")
    tibble::tibble(note = "CIF extraction failed; rerun with tidycmprsk >= 1.1.0")
  }
})

cat("Cumulative incidence curve data extracted:", nrow(cif_data), "rows.\n\n")

# --- 3.5 Analysis 3 forest plot figure ----------------------------------------
# Three models side by side for SAT and SBT primary exposure terms
# Discrete-time OR | Cox HR | Fine-Gray sHR

cat("-- 3.5 Analysis 3 forest plot figure\n")

# Assemble results for forest plot
forest_data_A3 <- bind_rows(
  tidy_dt_primary %>%
    filter(str_detect(term, "SAT_delivered_primary|SBT_delivered_2min")) %>%
    mutate(model = "Discrete-Time\n(Primary)", model_order = 1),
  tidy_cox %>%
    filter(str_detect(term, "SAT_delivered_primary|SBT_delivered_2min")) %>%
    mutate(model = "Cause-Specific\nCox", model_order = 2),
  tidy_fg %>%
    filter(str_detect(term, "SAT_prop_final_primary|SBT_prop_final_2min")) %>%
    mutate(
      model = "Fine-Gray\nSubdist.", model_order = 3,
      # Relabel Fine-Gray terms to match discrete-time for plotting
      term  = case_when(
        str_detect(term, "SAT") ~ "SAT_delivered_primary",
        str_detect(term, "SBT") ~ "SBT_delivered_2min",
        TRUE ~ term
      )
    )
) %>%
  mutate(
    trial       = if_else(str_detect(term, "SAT"), "SAT", "SBT"),
    model       = factor(model, levels = c("Discrete-Time\n(Primary)",
                                           "Cause-Specific\nCox",
                                           "Fine-Gray\nSubdist.")),
    est_label   = paste0(round(estimate, 2),
                         " (", round(conf.low, 2),
                         "–",  round(conf.high, 2), ")")
  )

fig_A3_forest <- ggplot(forest_data_A3,
                        aes(x = estimate, xmin = conf.low, xmax = conf.high,
                            y = model, color = trial)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_errorbarh(height = 0.2, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = est_label), hjust = -0.1, size = 3) +
  facet_wrap(~ trial, ncol = 1) +
  scale_color_manual(values = clr_trial) +
  scale_x_continuous(limits = c(0.5, 4.5)) +
  labs(
    title    = "Analysis 3: Time to Extubation -- Three-Model Comparison",
    subtitle = "OR (discrete-time) | HR (Cox) | sHR (Fine-Gray) with 95% CI",
    x        = "Effect Estimate (OR or HR; reference = 1)",
    y        = NULL,
    caption  = paste0(
      "Primary model: discrete-time logistic regression ",
      "(ns(vent_day, df=3); exact tied-event likelihood).\n",
      "Fine-Gray uses episode-level SAT/SBT proportion; ",
      "Cox and discrete-time use daily binary exposure.\n"
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_A3_forest, "A3_tte_outcomes/figures", "fig_A3_dt_forest.png",
           width = 9, height = 6)

# --- 3.S SENSITIVITY: Alternative exposure definitions -----------------------

cat("-- 3.S Sensitivity: Alternative exposure definitions\n")

run_dt_sensitivity <- function(sat_var, sbt_var, label, data) {
  cat("  Running sensitivity:", label, "\n")
  
  covariates_tv_s       <- drop_single_level(covariates_tv,       data)
  covariates_baseline_s <- drop_single_level(covariates_baseline, data)
  
  f_sens <- reformulate(
    termlabels = c(
      sat_var, sbt_var,
      "ns(vent_day, df = 3)",
      covariates_tv_s,
      covariates_baseline_s,
      if (!single_hospital) re_hosp
    ),
    response = "extubated"
  )
  
  sens_warnings <- character(0)
  
  fit_sens <- withCallingHandlers(
    if (single_hospital) {
      glm(f_sens, data = data, family = binomial(link = "logit"))
    } else {
      glmer(f_sens, data = data, family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl   = list(maxfun = 2e5)))
    },
    warning = function(w) {
      sens_warnings <<- c(sens_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  if (length(sens_warnings) > 0)
    cat("  WARNINGS:", paste(sens_warnings, collapse = "; "), "\n")
  
  tidy_out <- if (single_hospital) {
    broom::tidy(fit_sens, conf.int = TRUE, exponentiate = TRUE)
  } else {
    broom.mixed::tidy(fit_sens, conf.int = TRUE, exponentiate = TRUE,
                      effects = "fixed")
  }
  
  tidy_out %>%
    filter(str_detect(term, sat_var) | str_detect(term, sbt_var)) %>%
    mutate(sensitivity    = label,
           sat_exposure   = sat_var,
           sbt_exposure   = sbt_var)
}

sens_3S1 <- run_dt_sensitivity(
  "SAT_delivered_modified", "SBT_delivered_2min",
  "3S1_modified_SAT", df_dt
)
sens_3S2 <- run_dt_sensitivity(
  "SAT_delivered_primary", "SBT_delivered_5min",
  "3S2_5min_SBT", df_dt
)

results_3_sensitivity <- bind_rows(sens_3S1, sens_3S2)

cat("\nSensitivity results (Analysis 3):\n")
print(results_3_sensitivity %>%
        select(sensitivity, term, estimate, conf.low, conf.high, p.value))
cat("\n")

# --- 3.6 Analysis 3 exports ---------------------------------------------------

cat("-- 3.6 Exporting Analysis 3 outputs\n")

export_csv(
  tidy_dt_primary %>% mutate(model = "A3_primary_discrete_time",
                             warnings = paste(fit_dt_primary_warnings,
                                              collapse = "; ")),
  "A3_tte_outcomes/models", "A3_dt_primary_coefs.csv"
)

export_csv(
  re_var_dt,
  "A3_tte_outcomes/models", "A3_dt_primary_re_variance.csv"
)

export_csv(
  tidy_cox %>% mutate(model    = "A3_secondary_cox",
                      warnings = paste(fit_cox_warnings, collapse = "; ")),
  "A3_tte_outcomes/models", "A3_cox_secondary_coefs.csv"
)

export_csv(
  tidy_fg %>% mutate(model    = "A3_secondary_finegray",
                     warnings = paste(fit_fg_warnings, collapse = "; ")),
  "A3_tte_outcomes/models", "A3_fg_secondary_coefs.csv"
)

export_csv(
  cif_data,
  "A3_tte_outcomes/models", "A3_fg_cumulative_incidence.csv"
)

export_csv(
  cif_data,   # same data, figure-ready copy in figures/a3/
  "A3_tte_outcomes/figures", "A3_fig_cif_curves.csv"
)

export_csv(
  results_3_sensitivity %>% mutate(analysis = "A3"),
  "A3_tte_outcomes/models", "A3_sensitivity_coefs.csv"
)

export_rds(fit_dt_primary, "A3_tte_outcomes/models", "A3_fit_dt_primary.rds")
export_rds(fit_cox,        "A3_tte_outcomes/models", "A3_fit_cox_secondary.rds")
export_rds(fit_fg,         "A3_tte_outcomes/models", "A3_fit_fg_secondary.rds")

# --- 3.7 Patient-only GLM (no hospital RE) -- adjustment comparison input -----
# Fits the same discrete-time model without the hospital random intercept.
# Used exclusively by the Section 7 adjustment comparison figure.
# Exported so the CC can pool/audit both adjustment levels symmetrically.

cat("-- 3.7 Patient-only GLM (no hospital RE) for adjustment comparison\n")

f_dt_nore <- reformulate(
  termlabels = c(
    "SAT_delivered_primary", "SBT_delivered_2min",
    "ns(vent_day, df = 3)",
    covariates_tv_dt, covariates_baseline_dt
  ),
  response = "extubated"
)

fit_dt_nore_warnings <- character(0)

fit_dt_nore <- withCallingHandlers(
  glm(f_dt_nore, data = df_dt, family = binomial(link = "logit")),
  warning = function(w) {
    fit_dt_nore_warnings <<- c(fit_dt_nore_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_dt_nore_warnings) > 0)
  cat("WARNINGS:", paste(fit_dt_nore_warnings, collapse = "; "), "\n")

tidy_dt_nore <- broom::tidy(fit_dt_nore, conf.int = TRUE, exponentiate = TRUE)

cat("A3 patient-only GLM -- primary exposure terms:\n")
tidy_dt_nore %>%
  filter(str_detect(term, "SAT_delivered|SBT_delivered")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

export_csv(
  tidy_dt_nore %>% mutate(
    model      = "A3_dt_patient_only_glm",
    adjustment = "patient_only_no_hospital_RE",
    warnings   = paste(fit_dt_nore_warnings, collapse = "; ")
  ),
  "A3_tte_outcomes/models", "A3_dt_patient_only_coefs.csv"
)

cat("\nAnalysis 3 complete.\n\n")

# =============================================================================
# SECTION 4: ANALYSIS 4 -- VFD-28 TWO-PART MODEL
# =============================================================================

cat("============================================================\n")
cat("ANALYSIS 4: Ventilator-Free Days at 28 Days (VFD-28)\n")
cat("============================================================\n\n")

# NOTE: Two-part model replaces ZINB (per Hajage NEJM Evidence 2025)
# ZINB conflates mortality and liberation into one uninterpretable coefficient.
# Two-part separates them:
#   Part 1: Is SAT/SBT associated with being alive at 28 days?
#   Part 2: Among survivors, is SAT/SBT associated with more VFDs?
# Exposure: SAT_prop_final_primary, SBT_prop_final_2min (episode-level summary)
# Covariates: age, sex, CCI, SOFA_mean, FiO2_mean, PEEP_mean, sedation_mean
#        
#             NEE_mean excluded -- redundant with SOFA_mean cardiovascular component

# --- 4.0 Build Analysis 4 base dataset ----------------------------------------

df_a4 <- df_hosp %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean),                            # NEE_mean excluded -- redundant with SOFA_mean
    !is.na(alive_28d), !is.na(VFD_28)
  )

# Survivor subset for Part 2
df_a4_surv <- df_a4 %>% filter(survivor_28d == 1)

# Log complete case step to waterfall
waterfall <- log_step(waterfall, "5_complete_case_A4",
                      "File2", nrow(df_a4),
                      nrow(df_hosp) - nrow(df_a4),
                      "Complete case filter: Analysis 4 model covariates (NEE_mean excluded -- redundant with SOFA_mean)")

cat("Analysis 4 dataset:           ", nrow(df_a4), "hospitalizations\n")
cat("Alive at 28d (Part 1 outcome):", sum(df_a4$alive_28d), "\n")
cat("Survivors for Part 2:         ", nrow(df_a4_surv), "\n")
cat("Median VFD_28 [IQR]:          ",
    median(df_a4$VFD_28), "[",
    quantile(df_a4$VFD_28, 0.25), "-",
    quantile(df_a4$VFD_28, 0.75), "]\n\n")

# Shared covariate vector for Analysis 4 and 5
# hospital_type and location_type removed (Thread 4 amendment); NEE_mean removed (Q7)
covars_episode <- covariates_mean   # defined in Section 2.9

# --- 4.1 Part 1: Mixed-effects logistic -- alive at 28 days -------------------

cat("-- 4.1 Part 1: Alive at 28 days (mixed-effects logistic)\n")

# Apply single-level factor check to episode covariates
covars_episode_a4p1 <- drop_single_level(covars_episode, df_a4)

if (length(covars_episode_a4p1) < length(covars_episode)) {
  dropped_a4p1 <- setdiff(covars_episode, covars_episode_a4p1)
  cat("NOTE: Dropping single-level covariates from A4 Part 1 model:",
      paste(dropped_a4p1, collapse = ", "), "\n")
}

f_a4_part1 <- reformulate(
  termlabels = c(
    "SAT_prop_final_primary",
    "SBT_prop_final_2min",
    covars_episode_a4p1,
    if (!single_hospital) re_hosp
  ),
  response = "alive_28d"
)

cat("Formula:\n"); print(f_a4_part1); cat("\n")

fit_a4p1_warnings <- character(0)

fit_a4_part1 <- withCallingHandlers(
  if (single_hospital) {
    glm(f_a4_part1, data = df_a4, family = binomial(link = "logit"))
  } else {
    glmer(f_a4_part1, data = df_a4, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 2e5)))
  },
  warning = function(w) {
    fit_a4p1_warnings <<- c(fit_a4p1_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_a4p1_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_a4p1_warnings, collapse = "\n"), "\n")
}

tidy_a4p1 <- if (single_hospital) {
  broom::tidy(fit_a4_part1, conf.int = TRUE, exponentiate = TRUE)
} else {
  broom.mixed::tidy(fit_a4_part1, conf.int = TRUE, exponentiate = TRUE,
                    effects = "fixed")
}

cat("\nPart 1 results (OR for alive at 28d):\n")
tidy_a4p1 %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- 4.2 Part 2: Mixed-effects negative binomial -- VFDs among survivors ------
# Standard NB (not zero-truncated): survivors CAN have VFD = 0
# (ventilated all 28 days but alive -- this is not an impossible value)
# Survivor subgroup: death_flag == 0

cat("-- 4.2 Part 2: VFD-28 among survivors (mixed-effects NB)\n")

cat("Survivor dataset:             ", nrow(df_a4_surv), "hospitalizations\n")
cat("VFD_28 = 0 among survivors:  ",
    sum(df_a4_surv$VFD_28 == 0), "\n")
cat("Median VFD_28 [IQR] survivors:",
    median(df_a4_surv$VFD_28), "[",
    quantile(df_a4_surv$VFD_28, 0.25), "-",
    quantile(df_a4_surv$VFD_28, 0.75), "]\n\n")

# Apply single-level factor check to survivor subset
covars_episode_a4p2 <- drop_single_level(covars_episode, df_a4_surv)

if (length(covars_episode_a4p2) < length(covars_episode)) {
  dropped_a4p2 <- setdiff(covars_episode, covars_episode_a4p2)
  cat("NOTE: Dropping single-level covariates from A4 Part 2 model:",
      paste(dropped_a4p2, collapse = ", "), "\n")
}

f_a4_part2 <- reformulate(
  termlabels = c(
    "SAT_prop_final_primary",
    "SBT_prop_final_2min",
    covars_episode_a4p2,
    if (!single_hospital) re_hosp
  ),
  response = "VFD_28"
)

cat("Formula:\n"); print(f_a4_part2); cat("\n")

fit_a4p2_warnings <- character(0)

fit_a4_part2 <- withCallingHandlers(
  glmmTMB(f_a4_part2, data = df_a4_surv,
          family = nbinom2(link = "log")),
  warning = function(w) {
    fit_a4p2_warnings <<- c(fit_a4p2_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_a4p2_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_a4p2_warnings, collapse = "\n"), "\n")
}

tidy_a4p2 <- broom.mixed::tidy(fit_a4_part2, conf.int = TRUE,
                               exponentiate = TRUE, effects = "fixed")

cat("\nPart 2 results (IRR for VFD-28 among survivors):\n")
tidy_a4p2 %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- 4.3 Descriptive VFD-28 by delivery group --------------------------------
# Groups: neither / SAT-only / SBT-only / both
# Using primary exposure definitions
# prop > 0 = received at least one eligible day with delivery

cat("-- 4.3 Descriptive VFD-28 by delivery group\n")

df_a4_desc <- df_a4 %>%
  mutate(
    delivery_group = case_when(
      SAT_prop_final_primary > 0 & SBT_prop_final_2min > 0 ~ "Both",
      SAT_prop_final_primary > 0 & SBT_prop_final_2min == 0 ~ "SAT only",
      SAT_prop_final_primary == 0 & SBT_prop_final_2min > 0 ~ "SBT only",
      TRUE ~ "Neither"
    ),
    delivery_group = factor(delivery_group,
                            levels = c("Neither", "SAT only",
                                       "SBT only", "Both"))
  )

vfd_desc <- df_a4_desc %>%
  group_by(delivery_group) %>%
  summarise(
    n            = n(),
    median_vfd   = median(VFD_28, na.rm = TRUE),
    q1_vfd       = quantile(VFD_28, 0.25, na.rm = TRUE),
    q3_vfd       = quantile(VFD_28, 0.75, na.rm = TRUE),
    pct_vfd0     = round(mean(VFD_28 == 0,  na.rm = TRUE) * 100, 1),
    pct_vfd28    = round(mean(VFD_28 == 28, na.rm = TRUE) * 100, 1),
    pct_death    = round(mean(death_flag == 1, na.rm = TRUE) * 100, 1),
    .groups      = "drop"
  )

cat("\nVFD-28 descriptive by delivery group:\n")
print(vfd_desc)
cat("\n")

# --- 4.4 Analysis 4 two-part figure -------------------------------------------
# Two-panel: (A) Part 1 ORs alive at 28d; (B) Part 2 IRRs VFD survivors

cat("-- 4.4 Analysis 4 two-part figure\n")

# Assemble plot data for both parts
a4_plot_data <- bind_rows(
  tidy_a4p1 %>%
    filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(part = "Part 1: Alive at 28 Days\n(Mixed-Effects Logistic)",
           metric = "OR"),
  tidy_a4p2 %>%
    filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(part = "Part 2: VFD-28 Among Survivors\n(Mixed-Effects NB)",
           metric = "IRR")
) %>%
  mutate(
    trial     = if_else(str_detect(term, "SAT"), "SAT", "SBT"),
    est_label = paste0(round(estimate, 2),
                       " (", round(conf.low, 2),
                       "–",  round(conf.high, 2), ")")
  )

fig_A4 <- ggplot(a4_plot_data,
                 aes(x = estimate, xmin = conf.low, xmax = conf.high,
                     y = trial, color = trial)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_errorbarh(height = 0.25, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(label = est_label), hjust = -0.15, size = 3.2) +
  facet_wrap(~ part, ncol = 2, scales = "free_x") +
  scale_color_manual(values = clr_trial) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    title    = "Analysis 4: VFD-28 Two-Part Model Results",
    subtitle = "Part 1: OR for survival to 28d  |  Part 2: IRR for VFDs among survivors",
    x        = "Effect Estimate (reference = 1)",
    y        = NULL,
    caption  = paste0(
      "Two-part model replaces ZINB (Hajage NEJM Evidence 2025).\n",
      "SBT null in Part 2 reflects survivor bias -- not a model failure ",
      "(see Analysis 5 mortality results).\n"
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_A4, "A4_VFD_outcomes/figures", "fig_A4_twopart.png",
           width = 10, height = 5)

# VFD-28 distribution percentiles for CC figure construction
vfd_percentiles <- df_a4_desc %>%
  group_by(delivery_group) %>%
  reframe(
    percentile = 1:99,
    vfd_value  = quantile(VFD_28, probs = (1:99)/100, na.rm = TRUE)
  )

# --- 4.S SENSITIVITY: Alternative exposure definitions -----------------------

cat("-- 4.S Sensitivity: Alternative exposure definitions (Analysis 4)\n")

run_a4_sensitivity <- function(sat_var, sbt_var, label, data, data_surv) {
  cat("  Running sensitivity:", label, "\n")
  
  covars_ep_s1 <- drop_single_level(covars_episode, data)
  covars_ep_s2 <- drop_single_level(covars_episode, data_surv)
  
  f_p1 <- reformulate(
    termlabels = c(sat_var, sbt_var, covars_ep_s1,
                   if (!single_hospital) re_hosp),
    response = "alive_28d"
  )
  f_p2 <- reformulate(
    termlabels = c(sat_var, sbt_var, covars_ep_s2,
                   if (!single_hospital) re_hosp),
    response = "VFD_28"
  )
  
  w1 <- w2 <- character(0)
  
  fit_p1 <- withCallingHandlers(
    if (single_hospital) glm(f_p1, data = data, family = binomial())
    else glmer(f_p1, data = data, family = binomial(),
               control = glmerControl(optimizer = "bobyqa",
                                      optCtrl   = list(maxfun = 2e5))),
    warning = function(w) {
      w1 <<- c(w1, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  fit_p2 <- withCallingHandlers(
    glmmTMB(f_p2, data = data_surv, family = nbinom2(link = "log")),
    warning = function(w) {
      w2 <<- c(w2, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  t1 <- if (single_hospital) {
    broom::tidy(fit_p1, conf.int = TRUE, exponentiate = TRUE)
  } else {
    broom.mixed::tidy(fit_p1, conf.int = TRUE, exponentiate = TRUE,
                      effects = "fixed")
  }
  t2 <- broom.mixed::tidy(fit_p2, conf.int = TRUE, exponentiate = TRUE,
                          effects = "fixed")
  
  bind_rows(
    t1 %>%
      filter(str_detect(term, sat_var) | str_detect(term, sbt_var)) %>%
      mutate(part = "Part1_alive28d"),
    t2 %>%
      filter(str_detect(term, sat_var) | str_detect(term, sbt_var)) %>%
      mutate(part = "Part2_VFD_survivors")
  ) %>%
    mutate(sensitivity = label,
           warnings_p1 = paste(w1, collapse = "; "),
           warnings_p2 = paste(w2, collapse = "; "))
}

# Build sensitivity datasets using alternative prop variables
df_a4_sens <- df_hosp %>%
  filter(
    !is.na(SAT_prop_final_modified), !is.na(SBT_prop_final_5min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean), !is.na(alive_28d), !is.na(VFD_28)
  )

df_a4_surv_sens <- df_a4_sens %>% filter(survivor_28d == 1)

sens_4S1 <- run_a4_sensitivity(
  "SAT_prop_final_modified", "SBT_prop_final_2min",
  "4S1_modified_SAT", df_a4_sens, df_a4_surv_sens
)
sens_4S2 <- run_a4_sensitivity(
  "SAT_prop_final_primary", "SBT_prop_final_5min",
  "4S2_5min_SBT", df_a4_sens, df_a4_surv_sens
)

results_4_sensitivity <- bind_rows(sens_4S1, sens_4S2)

cat("\nSensitivity results (Analysis 4):\n")
print(results_4_sensitivity %>%
        select(sensitivity, part, term, estimate, conf.low, conf.high, p.value))
cat("\n")

# --- 4.5 Analysis 4 exports ---------------------------------------------------

cat("-- 4.5 Exporting Analysis 4 outputs\n")

export_csv(
  tidy_a4p1 %>% mutate(model    = "A4_part1_alive28d",
                       warnings = paste(fit_a4p1_warnings, collapse = "; ")),
  "A4_VFD_outcomes/models", "A4_part1_alive28d_coefs.csv"
)

export_csv(
  tidy_a4p2 %>% mutate(model    = "A4_part2_VFD_survivors",
                       warnings = paste(fit_a4p2_warnings, collapse = "; ")),
  "A4_VFD_outcomes/models", "A4_part2_vfd_survivors_coefs.csv"
)

export_csv(
  vfd_desc %>% mutate(model = "A4_descriptive_VFD28"),
  "A4_VFD_outcomes/models", "A4_vfd28_descriptive.csv"
)

export_csv(
  vfd_desc %>% mutate(model = "A4_descriptive_VFD28"),
  "A4_VFD_outcomes/tables", "vfd28_descriptive.csv"
)

export_csv(
  vfd_percentiles,
  "A4_VFD_outcomes/figures", "A4_fig_vfd_distribution.csv"
)

export_csv(
  results_4_sensitivity %>% mutate(analysis = "A4"),
  "A4_VFD_outcomes/models", "A4_sensitivity_coefs.csv"
)

export_rds(fit_a4_part1, "A4_VFD_outcomes/models", "A4_fit_part1_alive28d.rds")
export_rds(fit_a4_part2, "A4_VFD_outcomes/models", "A4_fit_part2_VFD_survivors.rds")

# --- 4.6 Patient-only GLMs (no hospital RE) -- adjustment comparison input ----

cat("-- 4.6 Patient-only GLMs (no hospital RE) for adjustment comparison\n")

# Part 1: binomial GLM
f_a4p1_nore <- reformulate(
  termlabels = c("SAT_prop_final_primary", "SBT_prop_final_2min",
                 covars_episode_a4p1),
  response = "alive_28d"
)

fit_a4p1_nore_warnings <- character(0)
fit_a4p1_nore <- withCallingHandlers(
  glm(f_a4p1_nore, data = df_a4, family = binomial(link = "logit")),
  warning = function(w) {
    fit_a4p1_nore_warnings <<- c(fit_a4p1_nore_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
if (length(fit_a4p1_nore_warnings) > 0)
  cat("  Part 1 WARNINGS:", paste(fit_a4p1_nore_warnings, collapse = "; "), "\n")

tidy_a4p1_nore <- broom::tidy(fit_a4p1_nore, conf.int = TRUE, exponentiate = TRUE)

# Part 2: NB GLM (standard NB, not zero-truncated -- GLM family limitation;
# noted in model_type column so CC is aware of the distributional difference)
f_a4p2_nore <- reformulate(
  termlabels = c("SAT_prop_final_primary", "SBT_prop_final_2min",
                 covars_episode_a4p2),
  response = "VFD_28"
)

fit_a4p2_nore_warnings <- character(0)
fit_a4p2_nore <- tryCatch(
  withCallingHandlers(
    MASS::glm.nb(f_a4p2_nore, data = df_a4_surv),
    warning = function(w) {
      fit_a4p2_nore_warnings <<- c(fit_a4p2_nore_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) {
    cat("  Part 2 MASS::glm.nb failed -- falling back to quasipoisson.",
        "\n  Error:", conditionMessage(e), "\n")
    fit_a4p2_nore_warnings <<- c(fit_a4p2_nore_warnings,
                                  paste("glm.nb failed; quasipoisson used:",
                                        conditionMessage(e)))
    glm(f_a4p2_nore, data = df_a4_surv, family = quasipoisson(link = "log"))
  }
)
if (length(fit_a4p2_nore_warnings) > 0)
  cat("  Part 2 WARNINGS:", paste(fit_a4p2_nore_warnings, collapse = "; "), "\n")

tidy_a4p2_nore <- broom::tidy(fit_a4p2_nore, conf.int = TRUE, exponentiate = TRUE)

cat("A4 patient-only GLMs -- primary exposure terms:\n")
bind_rows(
  tidy_a4p1_nore %>% filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(part = "Part1_alive28d"),
  tidy_a4p2_nore %>% filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(part = "Part2_VFD_survivors")
) %>% select(part, term, estimate, conf.low, conf.high, p.value) %>% print()
cat("\n")

export_csv(
  bind_rows(
    tidy_a4p1_nore %>%
      mutate(model      = "A4_part1_alive28d_patient_only_glm",
             part       = "Part1_alive28d",
             adjustment = "patient_only_no_hospital_RE",
             model_type = "binomial_glm",
             warnings   = paste(fit_a4p1_nore_warnings, collapse = "; ")),
    tidy_a4p2_nore %>%
      mutate(model      = "A4_part2_VFD_patient_only_glm",
             part       = "Part2_VFD_survivors",
             adjustment = "patient_only_no_hospital_RE",
             model_type = if (inherits(fit_a4p2_nore, "negbin"))
                            "negative_binomial_glm" else "quasipoisson_glm",
             warnings   = paste(fit_a4p2_nore_warnings, collapse = "; "))
  ),
  "A4_VFD_outcomes/models", "A4_patient_only_coefs.csv"
)

cat("\nAnalysis 4 complete.\n\n")
# =============================================================================
# SECTION 5: ANALYSIS 5 -- ICU LOS AND IN-HOSPITAL MORTALITY
# =============================================================================

cat("============================================================\n")
cat("ANALYSIS 5: ICU LOS and In-Hospital Mortality\n")
cat("============================================================\n\n")


# SBT is strongly protective for mortality (Analysis 5.2) but paradoxically
# associated with LONGER ICU LOS (IRR > 1). This is a structural finding,
# not a model error. SBT recipients survive longer → longer ICU stays by
# definition. Confirmed via noNEE sensitivity run (n=8,726)
#
# Covariates: age, sex, CCI, SOFA_mean, FiO2_mean, PEEP_mean, sedation_mean
#             hospital_type/location_type removed -- collinear with
#             (1|hospital_id)
#             NEE_mean excluded -- redundant with SOFA_mean cardiovascular component

# --- 5.0 Build Analysis 5 base dataset ----------------------------------------

df_a5 <- df_hosp %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean),                            # NEE_mean excluded -- see note above
    !is.na(ICU_LOS), !is.na(death_flag),
    ICU_LOS >= 1    # Structural constraint: LOS >= 1 by cohort definition
  )

# Coerce ICU_LOS to integer via ceiling() -- required by glmmTMB's
# truncated_nbinom2 family, which expects strict integer counts and warns
# otherwise. Ceiling is the most defensible clinical convention for partial-
# day LOS (any time in the ICU on a calendar day = 1 bed day), and preserves
# the >= 1 constraint for all rows passing the filter above.
# Scoped to df_a5 only: df_hosp retains fractional values for Table 1
# descriptives and other analyses that do not model LOS as a count outcome.
n_los_noninteger <- sum(df_a5$ICU_LOS != floor(df_a5$ICU_LOS), na.rm = TRUE)
if (n_los_noninteger > 0) {
  cat("ICU_LOS ceiling coercion:", n_los_noninteger, "of", nrow(df_a5),
      "rows had fractional LOS (",
      round(n_los_noninteger / nrow(df_a5) * 100, 1),
      "%). Applied as.integer(ceiling(ICU_LOS)) to resolve",
      "truncated_nbinom2 warning.\n\n")
} else {
  cat("ICU_LOS: all values already integer. Ceiling coercion applied",
      "(no-op -- warning would not have fired).\n\n")
}
df_a5 <- df_a5 %>% mutate(ICU_LOS = as.integer(ceiling(ICU_LOS)))
waterfall <- log_step(waterfall, "5_complete_case_A5",
                      "File2", nrow(df_a5),
                      nrow(df_hosp) - nrow(df_a5),
                      "Complete case filter: Analysis 5 model covariates (NEE_mean excluded -- redundant with SOFA_mean)")

cat("Analysis 5 dataset:       ", nrow(df_a5), "hospitalizations\n")
cat("Deaths:                   ", sum(df_a5$death_flag), "\n")
cat("Median ICU LOS [IQR]:     ",
    median(df_a5$ICU_LOS), "[",
    quantile(df_a5$ICU_LOS, 0.25), "-",
    quantile(df_a5$ICU_LOS, 0.75), "]\n\n")

# --- 5.1 ICU LOS: Zero-truncated negative binomial ---------------------------
# ZTNB conditions on LOS >= 1 (correct -- zero LOS impossible in this cohort)
# Standard NB assigns probability mass to LOS = 0 -- incorrect for this data
# hospital_type/location_type removed as covariates (Thread 4 amendment)
# NEE_mean removed (Q7)

cat("-- 5.1 ICU LOS: Zero-truncated negative binomial\n")

# Apply single-level factor check
covars_episode_a5los <- drop_single_level(covars_episode, df_a5)

if (length(covars_episode_a5los) < length(covars_episode)) {
  dropped_a5los <- setdiff(covars_episode, covars_episode_a5los)
  cat("NOTE: Dropping single-level covariates from A5 LOS model:",
      paste(dropped_a5los, collapse = ", "), "\n")
}

f_a5_los <- reformulate(
  termlabels = c(
    "SAT_prop_final_primary",
    "SBT_prop_final_2min",
    covars_episode_a5los,
    if (!single_hospital) re_hosp
  ),
  response = "ICU_LOS"
)

cat("Formula:\n"); print(f_a5_los); cat("\n")

fit_a5los_warnings <- character(0)

fit_a5_los <- withCallingHandlers(
  glmmTMB(f_a5_los, data = df_a5,
          family = truncated_nbinom2(link = "log")),
  warning = function(w) {
    fit_a5los_warnings <<- c(fit_a5los_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_a5los_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_a5los_warnings, collapse = "\n"), "\n")
}

tidy_a5_los <- broom.mixed::tidy(fit_a5_los, conf.int = TRUE,
                                 exponentiate = TRUE, effects = "fixed")

cat("\nICU LOS results (IRR):\n")
tidy_a5_los %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# Overdispersion check
# Threshold: > 1.5 flagged as structurally motivated (bimodal LOS trajectories)
# Per running edit list item 3 -- threshold lowered from > 2 to > 1.5
pearson_resid    <- residuals(fit_a5_los, type = "pearson")
dispersion_stat  <- sum(pearson_resid^2) / df.residual(fit_a5_los)
dispersion_flag  <- case_when(
  dispersion_stat > 2.0 ~ "HIGH -- investigate",
  dispersion_stat > 1.5 ~ "MODERATE -- structurally expected in ICU LOS",
  TRUE                  ~ "ok"
)

cat("Overdispersion stat (Pearson chi2/df):", round(dispersion_stat, 3),
    "--", dispersion_flag, "\n")
cat("Note: Moderate overdispersion (1.5-2.0) is structurally motivated in\n")
cat("  mechanically ventilated cohorts due to bimodal LOS trajectories\n")
cat("  (fast-extubation vs. prolonged ventilation). Pilot stat: 1.65.\n\n")

# Overdispersion as standalone export
overdispersion_out <- tibble(
  model              = "A5_ICU_LOS_ZTNB",
  dispersion_stat    = round(dispersion_stat, 4),
  n_obs              = nrow(df_a5),
  df_residual        = df.residual(fit_a5_los),
  flag               = dispersion_flag,
  icu_los_coercion   = "ceiling",
  n_los_ceiled       = n_los_noninteger,
  pct_los_ceiled     = round(n_los_noninteger / nrow(df_a5) * 100, 2),
  note               = paste0(
    "Threshold: >1.5 = structurally expected (bimodal ICU LOS); ",
    ">2.0 = investigate. ",
    "Structural overdispersion reflects mixture of fast-extubation ",
    "and prolonged ventilation trajectories. ",
    "ICU_LOS coerced to integer via ceiling() before model fit; ",
    n_los_noninteger, " rows (", round(n_los_noninteger/nrow(df_a5)*100,1),
    "%) had fractional values in the source data."
  )
)

# SBT paradox note for LOS
cat("NOTE -- SBT LOS paradox: IRR > 1 for SBT is a structural finding.\n")
cat("  SBT recipients survive longer → longer ICU stays by definition.\n")
cat("  Confirmed: noNEE sensitivity run showed IRR strengthened not attenuated.\n")
cat("  Acknowledged as survivor bias per Q8 team decision. No model change.\n\n")

# --- 5.2 Mortality: Mixed-effects logistic regression ------------------------
# Model-based SEs only -- no robust/sandwich SEs on top of random intercepts
# Double-counting per Cameron & Miller 2015 (see SAP Section 6)
# hospital_type/location_type removed as covariates (Thread 4 amendment)
# NEE_mean removed (Q7)

cat("-- 5.2 Mortality: Mixed-effects logistic regression\n")

# Apply single-level factor check
covars_episode_a5mort <- drop_single_level(covars_episode, df_a5)

if (length(covars_episode_a5mort) < length(covars_episode)) {
  dropped_a5mort <- setdiff(covars_episode, covars_episode_a5mort)
  cat("NOTE: Dropping single-level covariates from A5 mortality model:",
      paste(dropped_a5mort, collapse = ", "), "\n")
}

f_a5_mort <- reformulate(
  termlabels = c(
    "SAT_prop_final_primary",
    "SBT_prop_final_2min",
    covars_episode_a5mort,
    if (!single_hospital) re_hosp
  ),
  response = "death_flag"
)

cat("Formula:\n"); print(f_a5_mort); cat("\n")

fit_a5mort_warnings <- character(0)

fit_a5_mort <- withCallingHandlers(
  if (single_hospital) {
    glm(f_a5_mort, data = df_a5, family = binomial(link = "logit"))
  } else {
    glmer(f_a5_mort, data = df_a5, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 2e5)))
  },
  warning = function(w) {
    fit_a5mort_warnings <<- c(fit_a5mort_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_a5mort_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_a5mort_warnings, collapse = "\n"), "\n")
}

tidy_a5_mort <- if (single_hospital) {
  broom::tidy(fit_a5_mort, conf.int = TRUE, exponentiate = TRUE)
} else {
  broom.mixed::tidy(fit_a5_mort, conf.int = TRUE, exponentiate = TRUE,
                    effects = "fixed")
}

cat("\nMortality results (adjusted OR):\n")
tidy_a5_mort %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- 5.3 Analysis 5 two-panel figure ------------------------------------------
# Panel A: IRRs for ICU LOS from ZTNB
# Panel B: ORs for mortality from logistic
# SBT LOS paradox annotated in caption per Q8

cat("-- 5.3 Analysis 5 figure\n")

a5_plot_data <- bind_rows(
  tidy_a5_los %>%
    filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(outcome = "Panel A: ICU LOS\n(ZTNB -- IRR)",
           metric  = "IRR"),
  tidy_a5_mort %>%
    filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(outcome = "Panel B: In-Hospital Mortality\n(Mixed Logistic -- OR)",
           metric  = "OR")
) %>%
  mutate(
    trial     = if_else(str_detect(term, "SAT"), "SAT", "SBT"),
    est_label = paste0(round(estimate, 2),
                       " (", round(conf.low, 2),
                       "–",  round(conf.high, 2), ")")
  )

fig_A5 <- ggplot(a5_plot_data,
                 aes(x = estimate, xmin = conf.low, xmax = conf.high,
                     y = trial, color = trial)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_errorbarh(height = 0.25, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(label = est_label), hjust = -0.15, size = 3.2) +
  facet_wrap(~ outcome, ncol = 2, scales = "free_x") +
  scale_color_manual(values = clr_trial) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    title    = "Analysis 5: ICU LOS and In-Hospital Mortality",
    subtitle = "IRR for ICU LOS (ZTNB)  |  OR for mortality (mixed logistic)",
    x        = "Effect Estimate (reference = 1)",
    y        = NULL,
    caption  = paste0(
      "Panel A: SBT IRR > 1 reflects survivor bias -- SBT recipients ",
      "survive longer \u2192 longer ICU stays by definition.\n",
      "Confirmed: sensitivity run on larger population showed IRR ",
      "strengthened (1.29), ruling out population composition.\n",
      "Overdispersion stat: ", round(dispersion_stat, 2),
      " (", dispersion_flag, "). ",
      "Model-based SEs only; no sandwich SEs (Cameron & Miller 2015).\n"
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_A5, "A5_mort_outcomes/figures", "fig_A5_los_mortality.png",
           width = 10, height = 5)

# LOS distribution percentiles for CC figure construction
los_percentiles <- df_a5 %>%
  mutate(
    delivery_group = case_when(
      SAT_prop_final_primary > 0 & SBT_prop_final_2min > 0 ~ "Both",
      SAT_prop_final_primary > 0 & SBT_prop_final_2min == 0 ~ "SAT only",
      SAT_prop_final_primary == 0 & SBT_prop_final_2min > 0 ~ "SBT only",
      TRUE ~ "Neither"
    )
  ) %>%
  group_by(delivery_group) %>%
  reframe(
    percentile = 1:99,
    los_value  = quantile(ICU_LOS, probs = (1:99)/100, na.rm = TRUE)
  )

# --- 5.S SENSITIVITY: Alternative exposure definitions -----------------------

cat("-- 5.S Sensitivity: Alternative exposure definitions (Analysis 5)\n")

run_a5_sensitivity <- function(sat_var, sbt_var, label, data) {
  cat("  Running sensitivity:", label, "\n")
  
  covars_ep_s <- drop_single_level(covars_episode, data)
  
  f_los <- reformulate(
    termlabels = c(sat_var, sbt_var, covars_ep_s,
                   if (!single_hospital) re_hosp),
    response = "ICU_LOS"
  )
  f_mort <- reformulate(
    termlabels = c(sat_var, sbt_var, covars_ep_s,
                   if (!single_hospital) re_hosp),
    response = "death_flag"
  )
  
  wl <- wm <- character(0)
  
  fit_los <- withCallingHandlers(
    glmmTMB(f_los, data = data,
            family = truncated_nbinom2(link = "log")),
    warning = function(w) {
      wl <<- c(wl, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  fit_mort <- withCallingHandlers(
    if (single_hospital) {
      glm(f_mort, data = data, family = binomial())
    } else {
      glmer(f_mort, data = data, family = binomial(),
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl   = list(maxfun = 2e5)))
    },
    warning = function(w) {
      wm <<- c(wm, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  t_los <- broom.mixed::tidy(fit_los, conf.int = TRUE,
                             exponentiate = TRUE, effects = "fixed") %>%
    filter(str_detect(term, sat_var) | str_detect(term, sbt_var)) %>%
    mutate(outcome = "ICU_LOS", warnings = paste(wl, collapse = "; "))
  
  t_mort <- if (single_hospital) {
    broom::tidy(fit_mort, conf.int = TRUE, exponentiate = TRUE)
  } else {
    broom.mixed::tidy(fit_mort, conf.int = TRUE, exponentiate = TRUE,
                      effects = "fixed")
  }
  t_mort <- t_mort %>%
    filter(str_detect(term, sat_var) | str_detect(term, sbt_var)) %>%
    mutate(outcome = "mortality", warnings = paste(wm, collapse = "; "))
  
  bind_rows(t_los, t_mort) %>%
    mutate(sensitivity = label)
}

df_a5_sens <- df_a5 %>%
  filter(!is.na(SAT_prop_final_modified), !is.na(SBT_prop_final_5min))

sens_5S1 <- run_a5_sensitivity(
  "SAT_prop_final_modified", "SBT_prop_final_2min",
  "5S1_modified_SAT", df_a5_sens
)
sens_5S2 <- run_a5_sensitivity(
  "SAT_prop_final_primary", "SBT_prop_final_5min",
  "5S2_5min_SBT", df_a5_sens
)

results_5_sensitivity <- bind_rows(sens_5S1, sens_5S2)

cat("\nSensitivity results (Analysis 5):\n")
print(results_5_sensitivity %>%
        select(sensitivity, outcome, term, estimate,
               conf.low, conf.high, p.value))
cat("\n")

# --- 5.4 Analysis 5 exports ---------------------------------------------------

cat("-- 5.4 Exporting Analysis 5 outputs\n")

export_csv(
  tidy_a5_los %>% mutate(model    = "A5_ICU_LOS_ZTNB",
                         warnings = paste(fit_a5los_warnings, collapse = "; ")),
  "A5_mort_outcomes/models", "A5_icu_los_coefs.csv"
)

export_csv(
  overdispersion_out,
  "A5_mort_outcomes/models", "A5_icu_los_overdispersion.csv"
)

export_csv(
  tidy_a5_mort %>% mutate(model    = "A5_mortality_logistic",
                          warnings = paste(fit_a5mort_warnings, collapse = "; ")),
  "A5_mort_outcomes/models", "A5_mortality_coefs.csv"
)

export_csv(
  results_5_sensitivity %>% mutate(analysis = "A5"),
  "A5_mort_outcomes/models", "A5_sensitivity_coefs.csv"
)

export_csv(
  los_percentiles,
  "A5_mort_outcomes/figures", "A5_fig_los_distribution.csv"
)

export_rds(fit_a5_los,  "A5_mort_outcomes/models", "A5_fit_ICU_LOS_ZTNB.rds")
export_rds(fit_a5_mort, "A5_mort_outcomes/models", "A5_fit_mortality_logistic.rds")

# --- 5.5 Patient-only GLMs (no hospital RE) -- adjustment comparison input ----

cat("-- 5.5 Patient-only GLMs (no hospital RE) for adjustment comparison\n")

# ICU LOS: standard NB (not zero-truncated -- GLM family limitation;
# noted in model_type so CC is aware of distributional difference vs primary)
f_a5los_nore <- reformulate(
  termlabels = c("SAT_prop_final_primary", "SBT_prop_final_2min",
                 covars_episode_a5los),
  response = "ICU_LOS"
)

fit_a5los_nore_warnings <- character(0)
fit_a5los_nore <- tryCatch(
  withCallingHandlers(
    MASS::glm.nb(f_a5los_nore, data = df_a5),
    warning = function(w) {
      fit_a5los_nore_warnings <<- c(fit_a5los_nore_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) {
    cat("  LOS MASS::glm.nb failed -- falling back to quasipoisson.",
        "\n  Error:", conditionMessage(e), "\n")
    fit_a5los_nore_warnings <<- c(fit_a5los_nore_warnings,
                                   paste("glm.nb failed; quasipoisson used:",
                                         conditionMessage(e)))
    glm(f_a5los_nore, data = df_a5, family = quasipoisson(link = "log"))
  }
)
if (length(fit_a5los_nore_warnings) > 0)
  cat("  LOS WARNINGS:", paste(fit_a5los_nore_warnings, collapse = "; "), "\n")

tidy_a5los_nore <- broom::tidy(fit_a5los_nore, conf.int = TRUE,
                                exponentiate = TRUE)

# Mortality: binomial GLM
f_a5mort_nore <- reformulate(
  termlabels = c("SAT_prop_final_primary", "SBT_prop_final_2min",
                 covars_episode_a5mort),
  response = "death_flag"
)

fit_a5mort_nore_warnings <- character(0)
fit_a5mort_nore <- withCallingHandlers(
  glm(f_a5mort_nore, data = df_a5, family = binomial(link = "logit")),
  warning = function(w) {
    fit_a5mort_nore_warnings <<- c(fit_a5mort_nore_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
if (length(fit_a5mort_nore_warnings) > 0)
  cat("  Mortality WARNINGS:", paste(fit_a5mort_nore_warnings, collapse = "; "), "\n")

tidy_a5mort_nore <- broom::tidy(fit_a5mort_nore, conf.int = TRUE,
                                 exponentiate = TRUE)

cat("A5 patient-only GLMs -- primary exposure terms:\n")
bind_rows(
  tidy_a5los_nore  %>% filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(outcome = "ICU_LOS"),
  tidy_a5mort_nore %>% filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(outcome = "mortality")
) %>% select(outcome, term, estimate, conf.low, conf.high, p.value) %>% print()
cat("\n")

export_csv(
  bind_rows(
    tidy_a5los_nore %>%
      mutate(model      = "A5_ICU_LOS_patient_only_glm",
             outcome    = "ICU_LOS",
             adjustment = "patient_only_no_hospital_RE",
             model_type = if (inherits(fit_a5los_nore, "negbin"))
                            "negative_binomial_glm" else "quasipoisson_glm",
             warnings   = paste(fit_a5los_nore_warnings, collapse = "; ")),
    tidy_a5mort_nore %>%
      mutate(model      = "A5_mortality_patient_only_glm",
             outcome    = "mortality",
             adjustment = "patient_only_no_hospital_RE",
             model_type = "binomial_glm",
             warnings   = paste(fit_a5mort_nore_warnings, collapse = "; "))
  ),
  "A5_mort_outcomes/models", "A5_patient_only_coefs.csv"
)
# SECTION 6: SENSITIVITY ANALYSIS -- AGE < 65 SUBGROUP
# =============================================================================
# Repeats primary models from A3, A4, A5 restricted to patients aged < 65.
#
# Methodological note: age remains in all models even after restriction --
# it continues to vary within the under-65 subgroup and serves as a continuous
# confounder. This is a restricted-cohort SA, not a stratified analysis.
#
# Exposure and covariate definitions: identical to primary analyses.
# Outputs mirror primary analysis CSVs with "SA_age65_" prefix.
# =============================================================================

cat("============================================================\n")
cat("SENSITIVITY ANALYSIS: Age < 65 Subgroup\n")
cat("============================================================\n\n")

# --- Build age-restricted datasets -------------------------------------------
# age_u65 derived flag created in ABTRISE_01_setup.R Section 2.6

n_u65_hosp <- sum(df_hosp$age_u65 == 1, na.rm = TRUE)
n_u65_pp   <- sum(df_pp$age_u65   == 1, na.rm = TRUE)

cat("Age < 65 subgroup:\n")
cat("  Hospitalizations (File 2):", n_u65_hosp, "\n")
cat("  Person-days (File 1):     ", n_u65_pp, "\n\n")

if (n_u65_hosp < 50) {
  cat("WARNING: n < 50 in age < 65 subgroup. SA results may be unstable.\n")
  cat("  Proceeding but interpret with caution.\n\n")
}

# Hospitalization-level (A4, A5, Fine-Gray)
df_hosp_u65 <- df_hosp %>% filter(age_u65 == 1)

# Person-period (A3 discrete-time and Cox)
# Need last_day joined -- rebuild from df_pp_u65
df_pp_u65 <- df_pp %>% filter(age_u65 == 1)

last_day_u65 <- suppressWarnings(
  df_pp_u65 %>%
    group_by(hospitalization_id) %>%
    summarise(
      first_extub = min(vent_day[extubated == 1],  na.rm = TRUE),
      first_death = min(vent_day[died_today == 1], na.rm = TRUE),
      .groups = "drop"
    )
) %>%
  mutate(
    first_extub = if_else(is.infinite(first_extub), NA_real_, first_extub),
    first_death = if_else(is.infinite(first_death), NA_real_, first_death),
    last_day    = pmin(first_extub, first_death, 28, na.rm = TRUE)
  )

df_dt_u65 <- df_pp_u65 %>%
  left_join(last_day_u65 %>% select(hospitalization_id, last_day),
            by = "hospitalization_id") %>%
  filter(vent_day <= last_day) %>%
  filter(
    !is.na(SAT_delivered_primary), !is.na(SBT_delivered_2min),
    !is.na(SOFA_prior), !is.na(FiO2_prior), !is.na(PEEP_prior),
    !is.na(sedation_prior)
  )

cat("Age < 65 discrete-time dataset:", nrow(df_dt_u65), "person-day rows,",
    n_distinct(df_dt_u65$hospitalization_id), "hospitalizations\n\n")

# Episode-level complete-case datasets for A4/A5
df_a4_u65 <- df_hosp_u65 %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean), !is.na(alive_28d), !is.na(VFD_28)
  )

df_a4_surv_u65 <- df_a4_u65 %>% filter(survivor_28d == 1)

df_a5_u65 <- df_hosp_u65 %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean), !is.na(ICU_LOS), !is.na(death_flag),
    ICU_LOS >= 1
  ) %>%
  # Mirrors primary A5 coercion: ceiling() for truncated_nbinom2 compatibility
  mutate(ICU_LOS = as.integer(ceiling(ICU_LOS)))

cat("Age < 65 complete-case sizes:\n")
cat("  A4 dataset:       ", nrow(df_a4_u65), "| Survivors:", nrow(df_a4_surv_u65), "\n")
cat("  A5 dataset:       ", nrow(df_a5_u65), "\n\n")

# --- SA A3: Discrete-time logistic (age < 65) --------------------------------

cat("-- SA_age65 A3: Discrete-time logistic\n")

covariates_baseline_u65 <- drop_single_level(covariates_baseline, df_dt_u65)
covariates_tv_u65       <- drop_single_level(covariates_tv,       df_dt_u65)

f_dt_u65 <- reformulate(
  termlabels = c(
    "SAT_delivered_primary",
    "SBT_delivered_2min",
    "ns(vent_day, df = 3)",
    covariates_tv_u65,
    covariates_baseline_u65,
    if (!single_hospital) re_hosp
  ),
  response = "extubated"
)

sa_dt_u65_warnings <- character(0)

fit_sa_dt_u65 <- withCallingHandlers(
  if (single_hospital) {
    glm(f_dt_u65, data = df_dt_u65, family = binomial(link = "logit"))
  } else {
    glmer(f_dt_u65, data = df_dt_u65, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 2e5)))
  },
  warning = function(w) {
    sa_dt_u65_warnings <<- c(sa_dt_u65_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(sa_dt_u65_warnings) > 0) {
  cat("WARNINGS:", paste(sa_dt_u65_warnings, collapse = "; "), "\n")
}

tidy_sa_dt_u65 <- if (single_hospital) {
  broom::tidy(fit_sa_dt_u65, conf.int = TRUE, exponentiate = TRUE)
} else {
  broom.mixed::tidy(fit_sa_dt_u65, conf.int = TRUE, exponentiate = TRUE,
                    effects = "fixed")
}

cat("SA A3 discrete-time (age < 65) -- primary exposure terms:\n")
tidy_sa_dt_u65 %>%
  filter(str_detect(term, "SAT_delivered|SBT_delivered")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- SA A3 Cox: Cause-specific Cox (age < 65) --------------------------------
# Mirrors Section 3.3 -- same model specification, restricted to age < 65

cat("-- SA_age65 A3 Cox: Cause-specific Cox (counting process)\n")

df_cox_u65 <- df_dt_u65 %>%
  mutate(
    tstart   = vent_day - 1,
    tstop    = vent_day,
    event_cs = as.integer(extubated == 1 & died_today == 0)
  )

covariates_tv_cox_u65       <- drop_single_level(covariates_tv,       df_cox_u65)
covariates_baseline_cox_u65 <- drop_single_level(covariates_baseline, df_cox_u65)

if (length(covariates_baseline_cox_u65) < length(covariates_baseline)) {
  dropped_cox_u65 <- setdiff(covariates_baseline, covariates_baseline_cox_u65)
  cat("NOTE: Dropping single-level covariates from SA A3 Cox model:",
      paste(dropped_cox_u65, collapse = ", "), "\n")
}

f_cox_u65 <- reformulate(
  termlabels = c(
    "SAT_delivered_primary",
    "SBT_delivered_2min",
    covariates_tv_cox_u65,
    covariates_baseline_cox_u65,
    if (!single_hospital) "frailty(hospital_id, distribution = 'gamma')"
  ),
  response = "Surv(tstart, tstop, event_cs)"
)

cat("Formula:\n"); print(f_cox_u65); cat("\n")

fit_cox_u65_warnings <- character(0)

fit_cox_u65 <- withCallingHandlers(
  coxph(f_cox_u65, data = df_cox_u65, ties = "efron"),
  warning = function(w) {
    fit_cox_u65_warnings <<- c(fit_cox_u65_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_cox_u65_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_cox_u65_warnings, collapse = "\n"), "\n")
}

tidy_cox_u65 <- broom::tidy(fit_cox_u65, conf.int = TRUE, exponentiate = TRUE)

cat("\nSA A3 Cox (age < 65) -- primary exposure terms:\n")
tidy_cox_u65 %>%
  filter(str_detect(term, "SAT_delivered|SBT_delivered")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- SA A3 Fine-Gray: subdistribution hazard (age < 65) ----------------------
# Mirrors Section 3.4 -- same model specification, restricted to age < 65

cat("-- SA_age65 A3 Fine-Gray: subdistribution hazard\n")

df_fg_u65 <- df_hosp_u65 %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean)
    # NOTE: do NOT filter on !is.na(time_to_extubation) -- deaths and
    # censored patients have NA here and must be retained (see Section 3.4)
  ) %>%
  mutate(
    time_fg = case_when(
      extubation_flag == 1 ~ pmin(time_to_extubation, 28),
      death_flag == 1      ~ pmin(days_to_death, 28),
      TRUE                 ~ pmin(n_vent_days, 28)
    )
    # event_fg already constructed in Section 2.5 -- do not rebuild
  )

cat("SA A3 Fine-Gray (age < 65) dataset:", nrow(df_fg_u65), "patients\n")
cat("Event distribution:\n")
print(table(df_fg_u65$event_fg, useNA = "always"))
cat("\n")

covars_episode_fg_u65 <- drop_single_level(
  c("age", "sex", "CCI",
    "SOFA_mean", "FiO2_mean", "PEEP_mean", "sedation_mean"),
  df_fg_u65
)

if (length(covars_episode_fg_u65) < 7) {
  dropped_fg_u65 <- setdiff(
    c("age","sex","CCI","SOFA_mean","FiO2_mean","PEEP_mean","sedation_mean"),
    covars_episode_fg_u65
  )
  cat("NOTE: Dropping single-level covariates from SA A3 Fine-Gray model:",
      paste(dropped_fg_u65, collapse = ", "), "\n")
}

fg_rhs_u65 <- c(
  "SAT_prop_final_primary",
  "SBT_prop_final_2min",
  covars_episode_fg_u65,
  if (!single_hospital) "cluster(hospital_id)"
)

fg_formula_u65 <- as.formula(
  paste("Surv(time_fg, event_fg) ~", paste(fg_rhs_u65, collapse = " + "))
)

cat("Formula:\n"); print(fg_formula_u65); cat("\n")

fit_fg_u65_warnings <- character(0)

fit_fg_u65 <- withCallingHandlers(
  crr(fg_formula_u65, data = df_fg_u65, failcode = "extubated"),
  warning = function(w) {
    fit_fg_u65_warnings <<- c(fit_fg_u65_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(fit_fg_u65_warnings) > 0) {
  cat("WARNINGS (logged, model retained):\n")
  cat(paste("  -", fit_fg_u65_warnings, collapse = "\n"), "\n")
}

tidy_fg_u65 <- broom::tidy(fit_fg_u65, conf.int = TRUE, exponentiate = TRUE)

cat("SA A3 Fine-Gray (age < 65) -- primary exposure terms:\n")
tidy_fg_u65 %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# Cumulative incidence curve data (age < 65) -- mirrors Section 3.4
cif_data_u65 <- tryCatch({
  tidy_cif_u65 <- tidycmprsk::tidy(fit_fg_u65)
  as.data.frame(tidy_cif_u65)
}, error = function(e) {
  cat("NOTE: SA_age65 CIF extraction via tidy() failed:",
      conditionMessage(e), "\n")
  if (!is.null(fit_fg_u65$cmprsk)) {
    as.data.frame(fit_fg_u65$cmprsk$uftime)
  } else {
    tibble::tibble(note = "CIF extraction failed; rerun with tidycmprsk >= 1.1.0")
  }
})

cat("SA_age65 cumulative incidence curve data extracted:",
    nrow(cif_data_u65), "rows.\n\n")

# --- SA A4 Part 1: Alive at 28d (age < 65) -----------------------------------

cat("-- SA_age65 A4 Part 1: Alive at 28d\n")

covars_a4p1_u65 <- drop_single_level(covars_episode, df_a4_u65)

f_a4p1_u65 <- reformulate(
  termlabels = c("SAT_prop_final_primary", "SBT_prop_final_2min",
                 covars_a4p1_u65, if (!single_hospital) re_hosp),
  response = "alive_28d"
)

sa_a4p1_u65_warnings <- character(0)

fit_sa_a4p1_u65 <- withCallingHandlers(
  if (single_hospital) {
    glm(f_a4p1_u65, data = df_a4_u65, family = binomial(link = "logit"))
  } else {
    glmer(f_a4p1_u65, data = df_a4_u65, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 2e5)))
  },
  warning = function(w) {
    sa_a4p1_u65_warnings <<- c(sa_a4p1_u65_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

tidy_sa_a4p1_u65 <- if (single_hospital) {
  broom::tidy(fit_sa_a4p1_u65, conf.int = TRUE, exponentiate = TRUE)
} else {
  broom.mixed::tidy(fit_sa_a4p1_u65, conf.int = TRUE, exponentiate = TRUE,
                    effects = "fixed")
}

cat("SA A4 Part 1 (age < 65) -- primary exposure terms:\n")
tidy_sa_a4p1_u65 %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- SA A4 Part 2: VFD-28 among survivors (age < 65) ------------------------

cat("-- SA_age65 A4 Part 2: VFD-28 among survivors\n")
cat("  Survivors (age < 65):", nrow(df_a4_surv_u65), "\n")

covars_a4p2_u65 <- drop_single_level(covars_episode, df_a4_surv_u65)

f_a4p2_u65 <- reformulate(
  termlabels = c("SAT_prop_final_primary", "SBT_prop_final_2min",
                 covars_a4p2_u65, if (!single_hospital) re_hosp),
  response = "VFD_28"
)

sa_a4p2_u65_warnings <- character(0)

fit_sa_a4p2_u65 <- withCallingHandlers(
  glmmTMB(f_a4p2_u65, data = df_a4_surv_u65, family = nbinom2(link = "log")),
  warning = function(w) {
    sa_a4p2_u65_warnings <<- c(sa_a4p2_u65_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

tidy_sa_a4p2_u65 <- broom.mixed::tidy(fit_sa_a4p2_u65, conf.int = TRUE,
                                       exponentiate = TRUE, effects = "fixed")

cat("SA A4 Part 2 (age < 65) -- primary exposure terms:\n")
tidy_sa_a4p2_u65 %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- SA A5 LOS: ZTNB (age < 65) ----------------------------------------------

cat("-- SA_age65 A5 LOS: ZTNB\n")

covars_a5los_u65 <- drop_single_level(covars_episode, df_a5_u65)

f_a5los_u65 <- reformulate(
  termlabels = c("SAT_prop_final_primary", "SBT_prop_final_2min",
                 covars_a5los_u65, if (!single_hospital) re_hosp),
  response = "ICU_LOS"
)

sa_a5los_u65_warnings <- character(0)

fit_sa_a5los_u65 <- withCallingHandlers(
  glmmTMB(f_a5los_u65, data = df_a5_u65,
          family = truncated_nbinom2(link = "log")),
  warning = function(w) {
    sa_a5los_u65_warnings <<- c(sa_a5los_u65_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

tidy_sa_a5los_u65 <- broom.mixed::tidy(fit_sa_a5los_u65, conf.int = TRUE,
                                        exponentiate = TRUE, effects = "fixed")

cat("SA A5 LOS (age < 65) -- primary exposure terms:\n")
tidy_sa_a5los_u65 %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- SA A5 Mortality: mixed logistic (age < 65) ------------------------------

cat("-- SA_age65 A5 Mortality: mixed logistic\n")

covars_a5mort_u65 <- drop_single_level(covars_episode, df_a5_u65)

f_a5mort_u65 <- reformulate(
  termlabels = c("SAT_prop_final_primary", "SBT_prop_final_2min",
                 covars_a5mort_u65, if (!single_hospital) re_hosp),
  response = "death_flag"
)

sa_a5mort_u65_warnings <- character(0)

fit_sa_a5mort_u65 <- withCallingHandlers(
  if (single_hospital) {
    glm(f_a5mort_u65, data = df_a5_u65, family = binomial(link = "logit"))
  } else {
    glmer(f_a5mort_u65, data = df_a5_u65, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 2e5)))
  },
  warning = function(w) {
    sa_a5mort_u65_warnings <<- c(sa_a5mort_u65_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

tidy_sa_a5mort_u65 <- if (single_hospital) {
  broom::tidy(fit_sa_a5mort_u65, conf.int = TRUE, exponentiate = TRUE)
} else {
  broom.mixed::tidy(fit_sa_a5mort_u65, conf.int = TRUE, exponentiate = TRUE,
                    effects = "fixed")
}

cat("SA A5 Mortality (age < 65) -- primary exposure terms:\n")
tidy_sa_a5mort_u65 %>%
  filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()
cat("\n")

# --- SA_age65 figures ----------------------------------------------------------

cat("-- SA_age65: Building figures\n")

# -- SA A3 figure: three-model comparison, mirrors primary fig_A3_forest -----
# Discrete-time, Cox, and Fine-Gray for SAT/SBT primary exposure terms --
# same structure as the main-analysis Section 3.5 figure, age < 65 subgroup

sa_a3_plot_data <- bind_rows(
  tidy_sa_dt_u65 %>%
    filter(str_detect(term, "SAT_delivered_primary|SBT_delivered_2min")) %>%
    mutate(model = "Discrete-Time\n(Age <65)", model_order = 1),
  tidy_cox_u65 %>%
    filter(str_detect(term, "SAT_delivered_primary|SBT_delivered_2min")) %>%
    mutate(model = "Cause-Specific\nCox", model_order = 2),
  tidy_fg_u65 %>%
    filter(str_detect(term, "SAT_prop_final_primary|SBT_prop_final_2min")) %>%
    mutate(
      model = "Fine-Gray\nSubdist.", model_order = 3,
      # Relabel Fine-Gray terms to match discrete-time for plotting
      term  = case_when(
        str_detect(term, "SAT") ~ "SAT_delivered_primary",
        str_detect(term, "SBT") ~ "SBT_delivered_2min",
        TRUE ~ term
      )
    )
) %>%
  mutate(
    trial       = if_else(str_detect(term, "SAT"), "SAT", "SBT"),
    model       = factor(model, levels = c("Discrete-Time\n(Age <65)",
                                           "Cause-Specific\nCox",
                                           "Fine-Gray\nSubdist.")),
    est_label   = paste0(round(estimate, 2),
                         " (", round(conf.low, 2),
                         "–",  round(conf.high, 2), ")")
  )

fig_SA_A3 <- ggplot(sa_a3_plot_data,
                    aes(x = estimate, xmin = conf.low, xmax = conf.high,
                        y = model, color = trial)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_errorbarh(height = 0.2, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = est_label), hjust = -0.1, size = 3) +
  facet_wrap(~ trial, ncol = 1) +
  scale_color_manual(values = clr_trial) +
  scale_x_continuous(limits = c(0.5, 4.5)) +
  labs(
    title    = "Sensitivity Analysis (Age < 65): Time to Extubation -- Three-Model Comparison",
    subtitle = "OR (discrete-time) | HR (Cox) | sHR (Fine-Gray) with 95% CI -- age < 65 subgroup",
    x        = "Effect Estimate (OR or HR; reference = 1)",
    y        = NULL,
    caption  = paste0(
      "Age < 65 subgroup: ", n_u65_hosp, " hospitalizations, ",
      n_u65_pp, " person-days. Age retained as continuous covariate within ",
      "the restricted cohort (restricted-cohort SA, not a stratified ",
      "analysis). Mirrors primary Analysis 3 three-model specification ",
      "(Section 3.5).\n"
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_SA_A3, "A3_tte_outcomes/figures", "SA_age65_A3_forest.png",
           width = 9, height = 6)

# -- SA A4 figure: two-part model, mirrors primary fig_A4 structure -----------

sa_a4_plot_data <- bind_rows(
  tidy_sa_a4p1_u65 %>%
    filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(part = "Part 1: Alive at 28 Days\n(Mixed-Effects Logistic)",
           metric = "OR"),
  tidy_sa_a4p2_u65 %>%
    filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(part = "Part 2: VFD-28 Among Survivors\n(Mixed-Effects NB)",
           metric = "IRR")
) %>%
  mutate(
    trial     = if_else(str_detect(term, "SAT"), "SAT", "SBT"),
    est_label = paste0(round(estimate, 2),
                       " (", round(conf.low, 2),
                       "–",  round(conf.high, 2), ")")
  )

fig_SA_A4 <- ggplot(sa_a4_plot_data,
                    aes(x = estimate, xmin = conf.low, xmax = conf.high,
                        y = trial, color = trial)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_errorbarh(height = 0.25, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(label = est_label), hjust = -0.15, size = 3.2) +
  facet_wrap(~ part, ncol = 2, scales = "free_x") +
  scale_color_manual(values = clr_trial) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    title    = "Sensitivity Analysis (Age < 65): VFD-28 Two-Part Model",
    subtitle = "Part 1: OR for survival to 28d  |  Part 2: IRR for VFDs among survivors",
    x        = "Effect Estimate (reference = 1)",
    y        = NULL,
    caption  = paste0(
      "Age < 65 subgroup: ", nrow(df_a4_u65), " hospitalizations, ",
      nrow(df_a4_surv_u65), " survivors for Part 2. Mirrors primary ",
      "Analysis 4 model specification (Section 4.4); restricted-cohort SA, ",
      "not a stratified analysis.\n"
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_SA_A4, "A4_VFD_outcomes/figures", "SA_age65_A4_twopart.png",
           width = 10, height = 5)

# -- SA A5 figure: two-panel LOS/mortality, mirrors primary fig_A5 ------------
# Overdispersion check repeated here for the u65 LOS model (not computed
# elsewhere in the SA_age65 section above) so the figure caption isn't
# reporting an unchecked model assumption -- mirrors the Section 5.1 check.

pearson_resid_u65   <- residuals(fit_sa_a5los_u65, type = "pearson")
dispersion_stat_u65 <- sum(pearson_resid_u65^2) / df.residual(fit_sa_a5los_u65)
dispersion_flag_u65 <- case_when(
  dispersion_stat_u65 > 2.0 ~ "HIGH -- investigate",
  dispersion_stat_u65 > 1.5 ~ "MODERATE -- structurally expected in ICU LOS",
  TRUE                      ~ "ok"
)
cat("SA_age65 A5 LOS overdispersion stat (Pearson chi2/df):",
    round(dispersion_stat_u65, 3), "--", dispersion_flag_u65, "\n\n")

sa_a5_plot_data <- bind_rows(
  tidy_sa_a5los_u65 %>%
    filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(outcome = "Panel A: ICU LOS\n(ZTNB -- IRR)",
           metric  = "IRR"),
  tidy_sa_a5mort_u65 %>%
    filter(str_detect(term, "SAT_prop|SBT_prop")) %>%
    mutate(outcome = "Panel B: In-Hospital Mortality\n(Mixed Logistic -- OR)",
           metric  = "OR")
) %>%
  mutate(
    trial     = if_else(str_detect(term, "SAT"), "SAT", "SBT"),
    est_label = paste0(round(estimate, 2),
                       " (", round(conf.low, 2),
                       "–",  round(conf.high, 2), ")")
  )

fig_SA_A5 <- ggplot(sa_a5_plot_data,
                    aes(x = estimate, xmin = conf.low, xmax = conf.high,
                        y = trial, color = trial)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_errorbarh(height = 0.25, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(label = est_label), hjust = -0.15, size = 3.2) +
  facet_wrap(~ outcome, ncol = 2, scales = "free_x") +
  scale_color_manual(values = clr_trial) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    title    = "Sensitivity Analysis (Age < 65): ICU LOS and In-Hospital Mortality",
    subtitle = "IRR for ICU LOS (ZTNB)  |  OR for mortality (mixed logistic)",
    x        = "Effect Estimate (reference = 1)",
    y        = NULL,
    caption  = paste0(
      "Age < 65 subgroup: ", nrow(df_a5_u65), " hospitalizations. ",
      "Overdispersion stat: ", round(dispersion_stat_u65, 2),
      " (", dispersion_flag_u65, "). Mirrors primary Analysis 5 model ",
      "specification (Section 5.3); model-based SEs only ",
      "(Cameron & Miller 2015).\n"
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_SA_A5, "A5_mort_outcomes/figures", "SA_age65_A5_los_mortality.png",
           width = 10, height = 5)

cat("SA_age65 figures exported.\n\n")

# --- SA_age65 exports ---------------------------------------------------------

cat("-- SA_age65: Exporting outputs\n")

export_csv(
  tidy_sa_dt_u65 %>%
    mutate(model       = "SA_age65_A3_discrete_time",
           subgroup    = "age_lt65",
           n_subgroup  = n_distinct(df_dt_u65$hospitalization_id),
           warnings    = paste(sa_dt_u65_warnings, collapse = "; ")),
  "A3_tte_outcomes/models", "SA_age65_A3_dt_coefs.csv"
)

export_csv(
  tidy_cox_u65 %>%
    mutate(model       = "SA_age65_A3_cox",
           subgroup    = "age_lt65",
           n_subgroup  = n_distinct(df_cox_u65$hospitalization_id),
           warnings    = paste(fit_cox_u65_warnings, collapse = "; ")),
  "A3_tte_outcomes/models", "SA_age65_A3_cox_coefs.csv"
)

export_csv(
  tidy_fg_u65 %>%
    mutate(model       = "SA_age65_A3_finegray",
           subgroup    = "age_lt65",
           n_subgroup  = nrow(df_fg_u65),
           warnings    = paste(fit_fg_u65_warnings, collapse = "; ")),
  "A3_tte_outcomes/models", "SA_age65_A3_fg_coefs.csv"
)

export_csv(
  cif_data_u65 %>% mutate(subgroup = "age_lt65"),
  "A3_tte_outcomes/models", "SA_age65_A3_fg_cumulative_incidence.csv"
)

export_csv(
  cif_data_u65 %>% mutate(subgroup = "age_lt65"),   # figure-ready copy, mirrors Section 3.6
  "A3_tte_outcomes/figures", "SA_age65_A3_fig_cif_curves.csv"
)

export_csv(
  bind_rows(
    tidy_sa_a4p1_u65 %>%
      mutate(model = "SA_age65_A4_part1_alive28d", part = "Part1_alive28d",
             warnings = paste(sa_a4p1_u65_warnings, collapse = "; ")),
    tidy_sa_a4p2_u65 %>%
      mutate(model = "SA_age65_A4_part2_VFD_survivors", part = "Part2_VFD_survivors",
             warnings = paste(sa_a4p2_u65_warnings, collapse = "; "))
  ) %>% mutate(subgroup = "age_lt65", n_subgroup = nrow(df_a4_u65)),
  "A4_VFD_outcomes/models", "SA_age65_A4_coefs.csv"
)

export_csv(
  bind_rows(
    tidy_sa_a5los_u65 %>%
      mutate(model = "SA_age65_A5_ICU_LOS_ZTNB", outcome = "ICU_LOS",
             warnings = paste(sa_a5los_u65_warnings, collapse = "; ")),
    tidy_sa_a5mort_u65 %>%
      mutate(model = "SA_age65_A5_mortality_logistic", outcome = "mortality",
             warnings = paste(sa_a5mort_u65_warnings, collapse = "; "))
  ) %>% mutate(subgroup = "age_lt65", n_subgroup = nrow(df_a5_u65)),
  "A5_mort_outcomes/models", "SA_age65_A5_coefs.csv"
)

cat("SA_age65 outputs exported.\n\n")

# =============================================================================
# SECTION 7: ADJUSTMENT COMPARISON FIGURE
# Mirrors: Figure 2, Plotkin et al., JAMA Health Forum 2024
#   (DOI:10.1001/jamahealthforum.2024.0636)
#
# Assembles the adjustment comparison figure from patient-only GLM estimates
# already fitted and exported in sections 3.7, 4.6, and 5.5, alongside the
# primary GLMM estimates from sections 3.2, 4.1/4.2, and 5.1/5.2.
# No models are fitted here -- this section is figure assembly only.
#
# For A4 Part 2 and A5 LOS, the patient-only series uses standard NB
# (MASS::glm.nb) rather than zero-truncated NB (GLM family limitation);
# this is noted in the exported CSVs via the model_type column.
#
# For single-hospital sites, both series reflect GLM estimates (identical).
# =============================================================================

cat("============================================================\n")
cat("SECTION 7: Adjustment Comparison Figure (assembly only)\n")
cat("============================================================\n\n")

# ---------------------------------------------------------------------------
# 7.1 ASSEMBLE COMPARISON PLOT DATA
# ---------------------------------------------------------------------------
# Patient-only tidy objects: tidy_dt_nore, tidy_a4p1_nore, tidy_a4p2_nore,
#   tidy_a5los_nore, tidy_a5mort_nore  (fitted in sections 3.7, 4.6, 5.5)
# GLMM tidy objects: tidy_dt_primary, tidy_a4p1, tidy_a4p2,
#   tidy_a5_los, tidy_a5_mort  (fitted in sections 3.2, 4.1, 4.2, 5.1, 5.2)

cat("-- 7.1 Assembling comparison plot data\n\n")

outcome_levels <- c(
  "A3: Extubation\n(Disc-Time OR)",
  "A4: Alive 28d\n(OR)",
  "A4: VFDs\n(IRR*)",
  "A5: ICU LOS\n(IRR*)",
  "A5: Mortality\n(OR)"
)

pull_exposure <- function(tidy_df, outcome_label, adjust_label,
                          sat_pattern = "SAT", sbt_pattern = "SBT") {
  tidy_df %>%
    filter(str_detect(term, sat_pattern) | str_detect(term, sbt_pattern)) %>%
    mutate(
      outcome   = outcome_label,
      adjust    = adjust_label,
      trial     = if_else(str_detect(term, sat_pattern), "SAT", "SBT"),
      outcome_f = factor(outcome_label, levels = outcome_levels)
    ) %>%
    select(trial, outcome, outcome_f, adjust, estimate, conf.low, conf.high)
}

nore_rows <- bind_rows(
  pull_exposure(tidy_dt_nore,    outcome_levels[1], "Patient-adjusted",
                sat_pattern = "SAT_delivered", sbt_pattern = "SBT_delivered"),
  pull_exposure(tidy_a4p1_nore,  outcome_levels[2], "Patient-adjusted"),
  pull_exposure(tidy_a4p2_nore,  outcome_levels[3], "Patient-adjusted"),
  pull_exposure(tidy_a5los_nore, outcome_levels[4], "Patient-adjusted"),
  pull_exposure(tidy_a5mort_nore,outcome_levels[5], "Patient-adjusted")
)

glmm_rows <- bind_rows(
  pull_exposure(tidy_dt_primary, outcome_levels[1],
                "Patient- and hospital-adjusted",
                sat_pattern = "SAT_delivered", sbt_pattern = "SBT_delivered"),
  pull_exposure(tidy_a4p1,       outcome_levels[2],
                "Patient- and hospital-adjusted"),
  pull_exposure(tidy_a4p2,       outcome_levels[3],
                "Patient- and hospital-adjusted"),
  pull_exposure(tidy_a5_los,     outcome_levels[4],
                "Patient- and hospital-adjusted"),
  pull_exposure(tidy_a5_mort,    outcome_levels[5],
                "Patient- and hospital-adjusted")
)

adj_comp_data <- bind_rows(nore_rows, glmm_rows) %>%
  mutate(
    adjust_f = factor(adjust,
                      levels = c("Patient-adjusted",
                                 "Patient- and hospital-adjusted"))
  )

export_csv(adj_comp_data, "diagnostics", "fig_adj_comparison_data.csv")

# ---------------------------------------------------------------------------
# 7.2 ADJUSTMENT COMPARISON FIGURE
# ---------------------------------------------------------------------------

cat("-- 7.2 Building adjustment comparison figure\n\n")

clr_pt_only <- JAMA_COLORS[1]
clr_pt_hosp <- JAMA_COLORS[2]

adj_shapes <- c("Patient-adjusted"               = 16L,
                "Patient- and hospital-adjusted"  = 15L)
adj_colors <- c("Patient-adjusted"               = clr_pt_only,
                "Patient- and hospital-adjusted"  = clr_pt_hosp)

dodge_w <- 0.35

fig_adj_comparison <- ggplot(
    adj_comp_data,
    aes(x     = outcome_f,
        y     = estimate,
        ymin  = conf.low,
        ymax  = conf.high,
        color = adjust_f,
        shape = adjust_f,
        group = adjust_f)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray55", linewidth = 0.6) +
  geom_line(aes(group = adjust_f),
            position = position_dodge(width = dodge_w),
            linewidth = 0.55, alpha = 0.7) +
  geom_errorbar(width = 0.18, linewidth = 0.7,
                position = position_dodge(width = dodge_w)) +
  geom_point(size = 3.5,
             position = position_dodge(width = dodge_w)) +
  scale_color_manual(values = adj_colors, name = NULL) +
  scale_shape_manual(values = adj_shapes, name = NULL) +
  scale_y_log10(
    breaks = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0),
    labels = c("0.50", "0.75", "1.00", "1.25", "1.50", "2.00", "3.00")
  ) +
  facet_wrap(~ trial, ncol = 2, labeller = labeller(
    trial = c(SAT = "SAT (Spontaneous Awakening Trial)",
              SBT = "SBT (Spontaneous Breathing Trial)")
  )) +
  labs(
    title    = "Effect Estimates by Adjustment Level Across All Primary Outcomes",
    subtitle = "Patient-adjusted (GLM) vs. patient- and hospital-adjusted (GLMM with hospital random intercept)",
    x        = NULL,
    y        = "Effect Estimate (OR or IRR; log scale)",
    caption  = paste0(
      "* IRR for A4 VFDs and A5 ICU LOS. Patient-only series uses standard NB (MASS::glm.nb)",
      " as approximation for ZTNB; see A4_patient_only_coefs.csv and A5_patient_only_coefs.csv.\n",
      "Dashed line = null (reference = 1). Error bars = 95% CI (Wald for GLM; profile for GLMM).\n",
      if (single_hospital)
        "NOTE: Single-hospital site -- hospital RE not estimable; both series reflect GLM estimates (identical).\n"
      else
        paste0("Hospital RE included in patient+hospital model | n_hospitals = ", n_hospitals, ".\n"),
      "Reference: Plotkin et al., JAMA Health Forum 2024 (DOI:10.1001/jamahealthforum.2024.0636)."
    )
  ) +
  theme_abtrise() +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.9, "lines"),
    axis.text.x      = element_text(size = 8, lineheight = 0.9),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 10)
  )

export_png(fig_adj_comparison, "diagnostics", "fig_adj_comparison.png",
           width = 13, height = 6)

cat("  Adjustment comparison figure exported.\n\n")

# =============================================================================
# SESSION INFO
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "diagnostics", prefix_file("session_info_a345.txt")))
cat("Session info saved.\n")
cat("\n=== Analyses 3-4-5 script complete ===\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
