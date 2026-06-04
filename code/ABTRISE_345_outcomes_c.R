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
#   outputs/models/a3/    A3_dt_primary_coefs.csv, A3_dt_primary_re_variance.csv,
#                         A3_cox_secondary_coefs.csv, A3_fg_secondary_coefs.csv,
#                         A3_fg_cumulative_incidence.csv, A3_sensitivity_coefs.csv,
#                         SA_age65_A3_dt_coefs.csv,
#                         A3_fit_*.rds
#   outputs/models/a4/    A4_part1_alive28d_coefs.csv, A4_part2_vfd_survivors_coefs.csv,
#                         A4_vfd28_descriptive.csv, A4_sensitivity_coefs.csv,
#                         SA_age65_A4_coefs.csv,
#                         A4_fit_*.rds
#   outputs/models/a5/    A5_icu_los_coefs.csv, A5_icu_los_overdispersion.csv,
#                         A5_mortality_coefs.csv, A5_sensitivity_coefs.csv,
#                         SA_age65_A5_coefs.csv,
#                         A5_fit_*.rds
#   outputs/figures/a3/   fig_A3_dt_forest.png, A3_fig_cif_curves.csv
#   outputs/figures/a4/   fig_A4_twopart.png, A4_fig_vfd_distribution.csv
#   outputs/figures/a5/   fig_A5_los_mortality.png, A5_fig_los_distribution.csv
#   outputs/tables/       vfd28_descriptive.csv
#   outputs/diagnostics/  session_info_a345.txt
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
    !is.na(sedation_prior),
    !is.na(hospital_type), !is.na(location_type)   # Q4
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
#                  + baseline covariates + hospital_type + location_type
#                  + (1|hospital_id)
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
# hospital_type and location_type as covariates (Q4)
# NEE_prior included as binary time-varying covariate

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
# hospital_type and location_type as covariates (Q4)
# NEE_mean excluded -- redundant with SOFA_mean

cat("-- 3.4 Secondary: Fine-Gray subdistribution hazard\n")
df_fg <- df_hosp %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(hospital_type), !is.na(location_type),
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
  c("age", "sex", "CCI", "hospital_type", "location_type",
    "SOFA_mean", "FiO2_mean", "PEEP_mean", "sedation_mean"),
  df_fg
)

if (length(covars_episode_fg) < 9) {
  dropped_fg <- setdiff(
    c("age","sex","CCI","hospital_type","location_type",
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
# tidycmprsk::crr() does not store $cuminc -- use tidycmprsk::tidy() to extract
# the cumulative incidence function at each observed event time.
# The resulting data frame has columns: time, outcome, estimate, std.error,
# conf.low, conf.high -- one row per time point per competing event type.
# This is used for the CIF figure at the CC (delivered/not-delivered groups
# cannot be reconstructed here since crr() uses episode-level proportions,
# not a binary group variable; CC constructs group-stratified CIF separately).

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
      "Cox and discrete-time use daily binary exposure."
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_A3_forest, "figures/a3", "fig_A3_dt_forest.png",
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
  "models/a3", "A3_dt_primary_coefs.csv"
)

export_csv(
  re_var_dt,
  "models/a3", "A3_dt_primary_re_variance.csv"
)

export_csv(
  tidy_cox %>% mutate(model    = "A3_secondary_cox",
                      warnings = paste(fit_cox_warnings, collapse = "; ")),
  "models/a3", "A3_cox_secondary_coefs.csv"
)

export_csv(
  tidy_fg %>% mutate(model    = "A3_secondary_finegray",
                     warnings = paste(fit_fg_warnings, collapse = "; ")),
  "models/a3", "A3_fg_secondary_coefs.csv"
)

export_csv(
  cif_data,
  "models/a3", "A3_fg_cumulative_incidence.csv"
)

export_csv(
  cif_data,   # same data, figure-ready copy in figures/a3/
  "figures/a3", "A3_fig_cif_curves.csv"
)

export_csv(
  results_3_sensitivity %>% mutate(analysis = "A3"),
  "models/a3", "A3_sensitivity_coefs.csv"
)

export_rds(fit_dt_primary, "models/a3", "A3_fit_dt_primary.rds")
export_rds(fit_cox,        "models/a3", "A3_fit_cox_secondary.rds")
export_rds(fit_fg,         "models/a3", "A3_fit_fg_secondary.rds")

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
# This allows clinical readers to see whether the benefit is from reduced
# mortality, faster liberation, or both.
#
# Exposure: SAT_prop_final_primary, SBT_prop_final_2min (episode-level summary)
# Covariates: age, sex, CCI, hospital_type, location_type (Q4),
#             SOFA_mean, FiO2_mean, PEEP_mean, sedation_mean
#             NEE_mean excluded -- redundant with SOFA_mean cardiovascular component

# --- 4.0 Build Analysis 4 base dataset ----------------------------------------

df_a4 <- df_hosp %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(hospital_type), !is.na(location_type),    # Q4
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
# hospital_type and location_type included (Q4); NEE_mean removed (Q7)
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
      "(see Analysis 5 mortality results)."
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_A4, "figures/a4", "fig_A4_twopart.png",
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
    !is.na(hospital_type), !is.na(location_type),
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
  "models/a4", "A4_part1_alive28d_coefs.csv"
)

export_csv(
  tidy_a4p2 %>% mutate(model    = "A4_part2_VFD_survivors",
                       warnings = paste(fit_a4p2_warnings, collapse = "; ")),
  "models/a4", "A4_part2_vfd_survivors_coefs.csv"
)

export_csv(
  vfd_desc %>% mutate(model = "A4_descriptive_VFD28"),
  "models/a4", "A4_vfd28_descriptive.csv"
)

export_csv(
  vfd_desc %>% mutate(model = "A4_descriptive_VFD28"),
  "tables", "vfd28_descriptive.csv"
)

export_csv(
  vfd_percentiles,
  "figures/a4", "A4_fig_vfd_distribution.csv"
)

export_csv(
  results_4_sensitivity %>% mutate(analysis = "A4"),
  "models/a4", "A4_sensitivity_coefs.csv"
)

export_rds(fit_a4_part1, "models/a4", "A4_fit_part1_alive28d.rds")
export_rds(fit_a4_part2, "models/a4", "A4_fit_part2_VFD_survivors.rds")

cat("\nAnalysis 4 complete.\n\n")
# =============================================================================
# SECTION 5: ANALYSIS 5 -- ICU LOS AND IN-HOSPITAL MORTALITY
# =============================================================================

cat("============================================================\n")
cat("ANALYSIS 5: ICU LOS and In-Hospital Mortality\n")
cat("============================================================\n\n")

# NOTE -- SBT LOS PARADOX (Q8 team decision -- acknowledge as survivor bias):
# SBT is strongly protective for mortality (Analysis 5.2) but paradoxically
# associated with LONGER ICU LOS (IRR > 1). This is a structural finding,
# not a model error. SBT recipients survive longer → longer ICU stays by
# definition. Confirmed via noNEE sensitivity run (n=8,726): IRR strengthened
# to 1.29, ruling out population composition as the driver. Survivor bias
# explanation retained; annotated in fig_A5_los_mortality.png caption.
# No model change per Q8 team decision.
#
# NOTE -- OVERDISPERSION (structural):
# ICU LOS in mechanically ventilated patients reflects a mixture of clinical
# trajectories (fast-extubation vs. prolonged ventilation). Residual
# overdispersion after ZTNB is structurally motivated. Overdispersion stat
# flagged at threshold > 1.5 (not > 2). Pilot stat: 1.65.
#
# Covariates: age, sex, CCI, hospital_type, location_type (Q4),
#             SOFA_mean, FiO2_mean, PEEP_mean, sedation_mean
#             NEE_mean excluded -- redundant with SOFA_mean cardiovascular component

# --- 5.0 Build Analysis 5 base dataset ----------------------------------------

df_a5 <- df_hosp %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(hospital_type), !is.na(location_type),    # Q4
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean),                            # NEE_mean excluded -- see note above
    !is.na(ICU_LOS), !is.na(death_flag),
    ICU_LOS >= 1    # Structural constraint: LOS >= 1 by cohort definition
  )

# Log complete case step to waterfall
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
# hospital_type and location_type as covariates (Q4)
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
  model           = "A5_ICU_LOS_ZTNB",
  dispersion_stat = round(dispersion_stat, 4),
  n_obs           = nrow(df_a5),
  df_residual     = df.residual(fit_a5_los),
  flag            = dispersion_flag,
  note            = paste0(
    "Threshold: >1.5 = structurally expected (bimodal ICU LOS); ",
    ">2.0 = investigate. ",
    "Structural overdispersion reflects mixture of fast-extubation ",
    "and prolonged ventilation trajectories."
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
# hospital_type and location_type as covariates (Q4)
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
      "Model-based SEs only; no sandwich SEs (Cameron & Miller 2015)."
    )
  ) +
  theme_abtrise() +
  theme(legend.position = "none")

export_png(fig_A5, "figures/a5", "fig_A5_los_mortality.png",
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
  "models/a5", "A5_icu_los_coefs.csv"
)

export_csv(
  overdispersion_out,
  "models/a5", "A5_icu_los_overdispersion.csv"
)

export_csv(
  tidy_a5_mort %>% mutate(model    = "A5_mortality_logistic",
                          warnings = paste(fit_a5mort_warnings, collapse = "; ")),
  "models/a5", "A5_mortality_coefs.csv"
)

export_csv(
  results_5_sensitivity %>% mutate(analysis = "A5"),
  "models/a5", "A5_sensitivity_coefs.csv"
)

export_csv(
  los_percentiles,
  "figures/a5", "A5_fig_los_distribution.csv"
)

export_rds(fit_a5_los,  "models/a5", "A5_fit_ICU_LOS_ZTNB.rds")
export_rds(fit_a5_mort, "models/a5", "A5_fit_mortality_logistic.rds")


# =============================================================================
# SECTION 6: SENSITIVITY ANALYSIS -- AGE < 65 SUBGROUP
# =============================================================================
# Repeats primary models from A3, A4, A5 restricted to patients aged < 65.
# Rationale: assess whether SAT/SBT construct validity estimates are consistent
# in younger patients, who may differ in sedation burden, illness severity,
# and extubation trajectory from older ICU patients.
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
    !is.na(sedation_prior),
    !is.na(hospital_type), !is.na(location_type)
  )

cat("Age < 65 discrete-time dataset:", nrow(df_dt_u65), "person-day rows,",
    n_distinct(df_dt_u65$hospitalization_id), "hospitalizations\n\n")

# Episode-level complete-case datasets for A4/A5
df_a4_u65 <- df_hosp_u65 %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(hospital_type), !is.na(location_type),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean), !is.na(alive_28d), !is.na(VFD_28)
  )

df_a4_surv_u65 <- df_a4_u65 %>% filter(survivor_28d == 1)

df_a5_u65 <- df_hosp_u65 %>%
  filter(
    !is.na(SAT_prop_final_primary), !is.na(SBT_prop_final_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(hospital_type), !is.na(location_type),
    !is.na(SOFA_mean), !is.na(FiO2_mean), !is.na(PEEP_mean),
    !is.na(sedation_mean), !is.na(ICU_LOS), !is.na(death_flag),
    ICU_LOS >= 1
  )

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

# --- SA_age65 exports ---------------------------------------------------------

cat("-- SA_age65: Exporting outputs\n")

export_csv(
  tidy_sa_dt_u65 %>%
    mutate(model       = "SA_age65_A3_discrete_time",
           subgroup    = "age_lt65",
           n_subgroup  = n_distinct(df_dt_u65$hospitalization_id),
           warnings    = paste(sa_dt_u65_warnings, collapse = "; ")),
  "models/a3", "SA_age65_A3_dt_coefs.csv"
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
  "models/a4", "SA_age65_A4_coefs.csv"
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
  "models/a5", "SA_age65_A5_coefs.csv"
)

cat("SA_age65 outputs exported.\n\n")

# =============================================================================
# SESSION INFO
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "diagnostics", prefix_file("session_info_a345.txt")))
cat("Session info saved.\n")
cat("\n=== Analyses 3-4-5 script complete ===\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
