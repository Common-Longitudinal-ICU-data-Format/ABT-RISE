# =============================================================================
# ABT-RISE: Site-Level Analysis Script 4 of 4
# ANALYSIS 6 -- Hospital Benchmarking (Risk-Adjusted Delivery Rates)
#
# SCRIPTS IN THIS SERIES:
#   ABTRISE_01_setup.R          -- run directly to review data quality
#   ABTRISE_02_criterion.R
#   ABTRISE_345_outcomes.R
#   ABTRISE_06_benchmarking.R   <- YOU ARE HERE
#
# HOW TO RUN:
#   Open this file and click Source (or run source("ABTRISE_06_benchmarking.R"))
#   Setup runs automatically -- do NOT run ABTRISE_01_setup.R separately first.
#
# WHAT THIS SCRIPT PRODUCES:
#   outputs/models/a6/    A6_SAT_glmm_coefs.csv, A6_SBT_glmm_coefs.csv,
#                         A6_re_variance_mor.csv, A6_blups_risk_adj_rates.csv,
#                         A6_icc_patient_diagnostic.csv,
#                         A6_hospital_aggregate_summary.csv  [CC pooled model input]
#                         A6_SAT_ccc_results.csv, A6_SBT_ccc_results.csv [if flowsheet]
#                         A6_fit_*.rds
#   outputs/figures/a6/   fig_A6_SAT_caterpillar.png, fig_A6_SBT_caterpillar.png,
#                         fig_A6_SAT_funnel.png, fig_A6_SBT_funnel.png,
#                         A6_*_caterpillar_data.csv, A6_*_funnel_limits_data.csv,
#                         A6_*_bland_altman_data.csv [if flowsheet]
#   outputs/diagnostics/  session_info_a6.txt (updated exclusion waterfall)
#
# ANALYSIS 6 DESIGN NOTES:
#   - Two separate GLMMs: SAT delivery and SBT delivery (eligible vent-days)
#   - 2-level model: vent-days within hospitals (patient RE omitted for
#     consistency with A3-A5; ICC(patient) exported as diagnostic)
#   - Single-hospital fallback: same as A3-A5 (glm, no RE, NA placeholders)
#   - Dual output tracks:
#       Track 1 (local): GLMM coefficients + BLUPs + MOR/VPC
#       Track 2 (CC):    Hospital aggregate summary for pooled model fitting
#   - CCC + Bland-Altman gated on site_has_flowsheet_sat / _sbt
#     minimum 3 hospitals with flowsheet data required
#   - Caterpillar plot (main figure) + funnel plot (supplemental)
#
# COORDINATING CENTER: Rush
# =============================================================================

# --- Load setup (runs Sections 0-2 automatically) ----------------------------
source(here::here("code", "ABTRISE_01_setup_c.R"))

# SECTION 2-B: ANALYSIS 6 -- HOSPITAL BENCHMARKING
# =============================================================================

cat("============================================================\n")
cat("ANALYSIS 6: Hospital Benchmarking (Risk-Adjusted Delivery Rates)\n")
cat("============================================================\n\n")

# DESIGN NOTES (see script header for full rationale):
#   - Two separate GLMMs: SAT delivery | SBT delivery
#   - Unit of analysis: eligible vent-day (one row per eligible day)
#   - Outcome: delivered (0/1) among eligible days
#   - Covariates: age, sex, CCI, SOFA_prior, FiO2_prior, PEEP_prior,
#                 sedation_prior, NEE_prior (binary), hospital_type, location_type
#   - Random intercept: (1|hospital_id) -- suppressed for single-hospital sites
#   - 2-level model only (vent-days within hospitals); patient RE omitted
#     for consistency with A3–A5
#   - ICC(patient) computed as diagnostic from residuals -- always exported
#   - Dual output tracks:
#       Track 1 (local GLMM): coefficients + BLUPs + MOR/VPC -- multi-hospital only
#       Track 2 (CC aggregate): hospital summary CSV -- always exported
#   - CCC + Bland-Altman: gated on site_has_flowsheet_sat / _sbt
#     minimum 3 flowsheet hospitals required to run

# ---------------------------------------------------------------------------
# 6.0 BUILD ANALYSIS 6 DATASETS
# ---------------------------------------------------------------------------

cat("-- 6.0 Building Analysis 6 datasets\n\n")

# Covariate set for A6: same time-varying covariates as A3 discrete-time
# plus baseline covariates. Age/sex/CCI joined from File 2 in Section 2.4.
# hospital_type and location_type already factored in df_pp.

# Covariate set for A6: defined in ABTRISE_01_setup.R Section 2.9
# and available via source(). Reproduced here for reference only:
#   covariates_a6 <- c("age", "sex", "CCI",
#                      "SOFA_prior", "FiO2_prior", "PEEP_prior",
#                      "sedation_prior", "NEE_prior",
#                      "hospital_type", "location_type")

# --- A6 SAT dataset: eligible SAT days, complete case ----------------------
df_a6_sat <- df_pp %>%
  filter(SAT_eligible == 1L) %>%
  filter(
    !is.na(SAT_delivered_primary),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_prior), !is.na(FiO2_prior), !is.na(PEEP_prior),
    !is.na(sedation_prior),
    !is.na(hospital_type), !is.na(location_type)
    # NEE_prior: binary-recoded in Section 2.4 (NA=0) -- no NA filter needed
  ) %>%
  mutate(delivered = as.integer(SAT_delivered_primary))

# --- A6 SBT dataset: eligible SBT days, complete case ----------------------
df_a6_sbt <- df_pp %>%
  filter(SBT_eligible == 1L) %>%
  filter(
    !is.na(SBT_delivered_2min),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_prior), !is.na(FiO2_prior), !is.na(PEEP_prior),
    !is.na(sedation_prior),
    !is.na(hospital_type), !is.na(location_type)
  ) %>%
  mutate(delivered = as.integer(SBT_delivered_2min))

# Log to waterfall
waterfall <- log_step(waterfall, "6_complete_case_A6_SAT",
                      "File1", n_distinct(df_a6_sat$hospitalization_id),
                      n_distinct(df_pp %>% filter(SAT_eligible==1) %>%
                                   pull(hospitalization_id)) -
                        n_distinct(df_a6_sat$hospitalization_id),
                      "A6 SAT: eligible days, complete case on covariates")

waterfall <- log_step(waterfall, "6_complete_case_A6_SBT",
                      "File1", n_distinct(df_a6_sbt$hospitalization_id),
                      n_distinct(df_pp %>% filter(SBT_eligible==1) %>%
                                   pull(hospitalization_id)) -
                        n_distinct(df_a6_sbt$hospitalization_id),
                      "A6 SBT: eligible days, complete case on covariates")

cat("A6 SAT dataset:", nrow(df_a6_sat), "eligible days |",
    n_distinct(df_a6_sat$hospitalization_id), "hospitalizations |",
    n_distinct(df_a6_sat$hospital_id), "hospitals\n")
cat("  SAT delivery rate (raw):",
    round(mean(df_a6_sat$delivered) * 100, 1), "%\n\n")

cat("A6 SBT dataset:", nrow(df_a6_sbt), "eligible days |",
    n_distinct(df_a6_sbt$hospitalization_id), "hospitalizations |",
    n_distinct(df_a6_sbt$hospital_id), "hospitals\n")
cat("  SBT delivery rate (raw):",
    round(mean(df_a6_sbt$delivered) * 100, 1), "%\n\n")

# ---------------------------------------------------------------------------
# 6.1–6.2 GLMM HELPER: run one benchmarking model (SAT or SBT)
# ---------------------------------------------------------------------------
# Returns list with: tidy coefficients, RE variance/ICC, BLUPs,
#   hospital risk-adjusted rates, MOR, VPC, ICC(patient) diagnostic,
#   convergence warnings

run_a6_glmm <- function(trial_label, data) {

  cat("-- Fitting A6 GLMM:", trial_label, "\n")

  n_hosp_a6 <- n_distinct(data$hospital_id)

  # Apply single-level factor check (same helper as A3–A5)
  covars_a6_site <- drop_single_level(covariates_a6, data)

  if (length(covars_a6_site) < length(covariates_a6)) {
    dropped <- setdiff(covariates_a6, covars_a6_site)
    cat("  NOTE: Dropping single-level covariates:",
        paste(dropped, collapse = ", "), "\n")
  }

  f_a6 <- reformulate(
    termlabels = c(covars_a6_site,
                   if (!single_hospital) re_hosp),
    response = "delivered"
  )

  cat("  Formula:\n"); print(f_a6); cat("\n")

  fit_warnings <- character(0)

  fit <- withCallingHandlers(
    if (single_hospital) {
      glm(f_a6, data = data, family = binomial(link = "logit"))
    } else {
      glmer(f_a6, data = data, family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl   = list(maxfun = 2e5)))
    },
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (length(fit_warnings) > 0) {
    cat("  WARNINGS (logged, model retained):\n")
    cat(paste("    -", fit_warnings, collapse = "\n"), "\n")
  }

  # --- Fixed effects ---------------------------------------------------------
  tidy_out <- if (single_hospital) {
    broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  } else {
    broom.mixed::tidy(fit, conf.int = TRUE, exponentiate = TRUE,
                      effects = "fixed")
  }

  # --- RE variance, VPC/ICC, MOR --------------------------------------------
  if (!single_hospital) {
    vc        <- as.data.frame(VarCorr(fit))
    sigma2_u  <- vc$vcov[1]          # hospital-level variance
    sigma2_e  <- pi^2 / 3            # logistic residual variance (standard)
    vpc       <- sigma2_u / (sigma2_u + sigma2_e)
    mor       <- exp(sqrt(2 * sigma2_u) * 0.6745)

    re_stats <- tibble(
      trial          = trial_label,
      model_type     = "glmer",
      sigma2_hosp    = round(sigma2_u, 6),
      sd_hosp        = round(sqrt(sigma2_u), 6),
      vpc_icc        = round(vpc, 4),
      MOR            = round(mor, 4),
      n_hospitals    = n_hosp_a6,
      single_hosp    = FALSE,
      warnings       = paste(fit_warnings, collapse = "; ")
    )

    cat("  Hospital-level variance (sigma2_u):", round(sigma2_u, 4), "\n")
    cat("  VPC/ICC:", round(vpc, 4),
        "| MOR:", round(mor, 4), "\n\n")

  } else {
    re_stats <- tibble(
      trial          = trial_label,
      model_type     = "glm_fallback",
      sigma2_hosp    = NA_real_,
      sd_hosp        = NA_real_,
      vpc_icc        = NA_real_,
      MOR            = NA_real_,
      n_hospitals    = n_hosp_a6,
      single_hosp    = TRUE,
      warnings       = paste(fit_warnings, collapse = "; ")
    )
    cat("  Single-hospital site: RE variance = NA, MOR = NA\n\n")
  }

  # --- BLUPs and risk-adjusted hospital rates --------------------------------
  # BLUPs = empirical Bayes shrinkage estimates (ranef())
  # Risk-adjusted rate = plogis(intercept + BLUP_j) -- the predicted
  # probability for a "reference patient" at hospital j
  # (Rogers et al. 2013: standard for provider profiling)

  if (!single_hospital) {
    blups <- ranef(fit)$hospital_id[, 1, drop = TRUE]
    hosp_ids_model <- rownames(ranef(fit)$hospital_id)
    intercept <- fixef(fit)[["(Intercept)"]]

    blup_df <- tibble(
      trial             = trial_label,
      hospital_id       = hosp_ids_model,
      blup              = blups,
      blup_se           = sqrt(attr(ranef(fit, condVar = TRUE)$hospital_id,
                                    "postVar")[1, 1, ]),
      risk_adj_rate     = plogis(intercept + blups),
      risk_adj_rate_lo  = plogis(intercept + blups -
                                   1.96 * sqrt(attr(
                                     ranef(fit, condVar = TRUE)$hospital_id,
                                     "postVar")[1, 1, ])),
      risk_adj_rate_hi  = plogis(intercept + blups +
                                   1.96 * sqrt(attr(
                                     ranef(fit, condVar = TRUE)$hospital_id,
                                     "postVar")[1, 1, ]))
    ) %>%
      # Join raw delivery rate for comparison
      left_join(
        data %>%
          group_by(hospital_id) %>%
          summarise(n_eligible   = n(),
                    n_delivered  = sum(delivered),
                    raw_rate     = mean(delivered),
                    .groups = "drop"),
        by = "hospital_id"
      ) %>%
      arrange(risk_adj_rate)

    cat("  BLUPs computed for", nrow(blup_df), "hospitals.\n")
    cat("  Risk-adjusted rate range: [",
        round(min(blup_df$risk_adj_rate), 3), ",",
        round(max(blup_df$risk_adj_rate), 3), "]\n\n")

  } else {
    # Single-hospital: one row with NA BLUPs; raw rate still informative
    blup_df <- data %>%
      group_by(hospital_id) %>%
      summarise(n_eligible  = n(),
                n_delivered = sum(delivered),
                raw_rate    = mean(delivered),
                .groups     = "drop") %>%
      mutate(trial            = trial_label,
             blup             = NA_real_,
             blup_se          = NA_real_,
             risk_adj_rate    = NA_real_,
             risk_adj_rate_lo = NA_real_,
             risk_adj_rate_hi = NA_real_)

    cat("  Single-hospital site: BLUPs = NA\n\n")
  }

  # --- ICC(patient) diagnostic -----------------------------------------------
  # Computed from model residuals using the latent variable approach:
  # ICC(patient) approx = variance explained by patient / total variance
  # Estimated via pearson residual variance partitioned by patient
  # This is a diagnostic only -- no modeling decision made at site level
  # CC reviews across sites to assess whether 3-level model is warranted

  cat("  Computing ICC(patient) diagnostic...\n")

  resid_df <- data %>%
    mutate(.fitted  = predict(fit, type = "response"),
           .resid_p = (delivered - .fitted) /
                        sqrt(.fitted * (1 - .fitted) + 1e-8)) %>%
    group_by(hospitalization_id) %>%
    summarise(patient_mean_resid = mean(.resid_p, na.rm = TRUE),
              .groups = "drop")

  var_patient   <- var(resid_df$patient_mean_resid, na.rm = TRUE)
  var_total_approx <- var(
    data %>%
      mutate(.fitted = predict(fit, type = "response"),
             .resid_p = (delivered - .fitted) /
                          sqrt(.fitted * (1 - .fitted) + 1e-8)) %>%
      pull(.resid_p),
    na.rm = TRUE
  )

  icc_patient <- if (var_total_approx > 0)
    round(var_patient / var_total_approx, 4) else NA_real_

  cat("  ICC(patient) [diagnostic]:", icc_patient, "\n")
  if (!is.na(icc_patient) && icc_patient > 0.05)
    cat("  NOTE: ICC(patient) > 0.05 -- flag for CC review.",
        "3-level model may be warranted at pooled analysis.\n")
  cat("\n")

  icc_patient_out <- tibble(
    trial               = trial_label,
    icc_patient         = icc_patient,
    var_patient_resid   = round(var_patient, 6),
    var_total_resid     = round(var_total_approx, 6),
    note = if (!is.na(icc_patient) && icc_patient > 0.05)
      "ICC(patient) > 0.05 -- flag for CC review; 3-level model may be warranted"
    else if (!is.na(icc_patient))
      "ICC(patient) <= 0.05 -- 2-level model adequate"
    else
      "ICC(patient) could not be computed"
  )

  list(
    tidy         = tidy_out %>% mutate(trial = trial_label,
                                       warnings = paste(fit_warnings,
                                                        collapse = "; ")),
    re_stats     = re_stats,
    blups        = blup_df,
    icc_patient  = icc_patient_out,
    fit          = fit,
    fit_warnings = fit_warnings
  )
}

# ---------------------------------------------------------------------------
# 6.1 SAT BENCHMARKING GLMM
# ---------------------------------------------------------------------------

cat("============================================================\n")
cat("6.1 SAT Benchmarking GLMM\n")
cat("============================================================\n\n")

a6_sat <- run_a6_glmm("SAT", df_a6_sat)

# ---------------------------------------------------------------------------
# 6.2 SBT BENCHMARKING GLMM
# ---------------------------------------------------------------------------

cat("============================================================\n")
cat("6.2 SBT Benchmarking GLMM\n")
cat("============================================================\n\n")

a6_sbt <- run_a6_glmm("SBT", df_a6_sbt)

# ---------------------------------------------------------------------------
# 6.3 HOSPITAL AGGREGATE SUMMARY (CC-LEVEL POOLED MODEL INPUT)
# ---------------------------------------------------------------------------
# Track 2 output -- always exported regardless of single vs. multi-hospital
# One row per hospital per trial (long format)
# Covariate profile averaged across eligible vent-days for that hospital
# CC uses this to fit a pooled GLMM across all sites without individual data

cat("-- 6.3 Hospital aggregate summary (CC pooled model input)\n\n")

build_hosp_summary <- function(trial_label, data, raw_rate_col = "delivered") {

  data %>%
    group_by(hospital_id) %>%
    summarise(
      trial              = trial_label,
      site_id            = site_id,
      n_eligible_days    = n(),
      n_delivered        = sum(.data[[raw_rate_col]]),
      raw_rate           = round(mean(.data[[raw_rate_col]]), 4),
      n_hospitalizations = n_distinct(hospitalization_id),
      # Covariate profile: means across eligible vent-days
      age_mean           = round(mean(age,            na.rm = TRUE), 2),
      sex_pct_male       = round(mean(sex == "M",     na.rm = TRUE) * 100, 1),
      CCI_mean           = round(mean(CCI,            na.rm = TRUE), 2),
      SOFA_prior_mean    = round(mean(SOFA_prior,     na.rm = TRUE), 2),
      FiO2_prior_mean    = round(mean(FiO2_prior,     na.rm = TRUE), 4),
      PEEP_prior_mean    = round(mean(PEEP_prior,     na.rm = TRUE), 2),
      sedation_prior_pct = round(mean(sedation_prior, na.rm = TRUE) * 100, 1),
      NEE_prior_pct      = round(mean(NEE_prior,      na.rm = TRUE) * 100, 1),
      # Fixed hospital-level characteristics (should be constant within hospital)
      hospital_type      = first(as.character(hospital_type)),
      location_type      = first(as.character(location_type)),
      .groups            = "drop"
    ) %>%
    mutate(single_hospital_site = single_hospital)
}

hosp_summary_sat <- build_hosp_summary("SAT", df_a6_sat)
hosp_summary_sbt <- build_hosp_summary("SBT", df_a6_sbt)
hosp_summary_all <- bind_rows(hosp_summary_sat, hosp_summary_sbt)

cat("Hospital aggregate summary (SAT):\n")
print(hosp_summary_sat %>% select(hospital_id, n_eligible_days,
                                   n_delivered, raw_rate, n_hospitalizations))
cat("\nHospital aggregate summary (SBT):\n")
print(hosp_summary_sbt %>% select(hospital_id, n_eligible_days,
                                   n_delivered, raw_rate, n_hospitalizations))
cat("\n")

# ---------------------------------------------------------------------------
# 6.4 CCC + BLAND-ALTMAN: EHR VS. FLOWSHEET HOSPITAL-LEVEL AGREEMENT
# ---------------------------------------------------------------------------
# Runs only at sites with flowsheet data AND >= 3 hospitals with flowsheet
# Compares algorithm delivery RATE per hospital vs. flowsheet delivery RATE
# This is the hospital-level agreement analysis (distinct from day-level A2)
# CCC: Lin's concordance correlation -- captures correlation + absolute agreement
# Bland-Altman: reveals systematic bias (algorithm consistently over/underestimates)

cat("-- 6.4 CCC + Bland-Altman: EHR vs. flowsheet agreement\n\n")

run_a6_ccc <- function(trial_label, alg_var, flowsheet_var,
                        elig_var, data = df_pp) {

  # Hospital-level rates: algorithm vs. flowsheet
  # Only hospitals with >= 1 eligible day with non-missing flowsheet data
  hosp_rates <- data %>%
    filter(.data[[elig_var]] == 1L,
           !is.na(.data[[flowsheet_var]]),
           !is.na(.data[[alg_var]])) %>%
    group_by(hospital_id) %>%
    summarise(
      n_days        = n(),
      rate_alg      = mean(as.integer(.data[[alg_var]]),      na.rm = TRUE),
      rate_flowsheet = mean(as.integer(.data[[flowsheet_var]]), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_days >= 5)   # need at least 5 eligible days per hospital to
                          # compute a stable rate

  n_hosp_fs <- nrow(hosp_rates)
  cat("  ", trial_label, "-- hospitals with flowsheet data (>=5 days):",
      n_hosp_fs, "\n")

  if (n_hosp_fs < 3) {
    cat("  SKIP: fewer than 3 hospitals with flowsheet data.\n")
    cat("  Bland-Altman and CCC uninformative below this threshold.\n\n")
    return(NULL)
  }

  # --- Lin's CCC via epiR::epi.ccc() ----------------------------------------
  ccc_result <- epiR::epi.ccc(
    x = hosp_rates$rate_alg,
    y = hosp_rates$rate_flowsheet
  )

  ccc_out <- tibble(
    trial          = trial_label,
    alg_var        = alg_var,
    flowsheet_var  = flowsheet_var,
    n_hospitals    = n_hosp_fs,
    CCC            = round(ccc_result$rho.c$est,  4),
    CCC_lo         = round(ccc_result$rho.c$lower, 4),
    CCC_hi         = round(ccc_result$rho.c$upper, 4),
    pearson_r      = round(ccc_result$r,           4),
    cb_coef        = round(ccc_result$C.b,         4),
    mean_diff      = round(mean(hosp_rates$rate_alg -
                                  hosp_rates$rate_flowsheet), 4),
    sd_diff        = round(sd(hosp_rates$rate_alg -
                                hosp_rates$rate_flowsheet), 4),
    loa_lo         = round(mean(hosp_rates$rate_alg -
                                  hosp_rates$rate_flowsheet) -
                             1.96 * sd(hosp_rates$rate_alg -
                                         hosp_rates$rate_flowsheet), 4),
    loa_hi         = round(mean(hosp_rates$rate_alg -
                                  hosp_rates$rate_flowsheet) +
                             1.96 * sd(hosp_rates$rate_alg -
                                         hosp_rates$rate_flowsheet), 4)
  )

  cat("  CCC:", round(ccc_result$rho.c$est, 4),
      "(", round(ccc_result$rho.c$lower, 4), "-",
      round(ccc_result$rho.c$upper, 4), ")\n")
  cat("  Mean difference (alg - flowsheet):",
      round(ccc_out$mean_diff, 4),
      "| LOA: [", round(ccc_out$loa_lo, 4), ",",
      round(ccc_out$loa_hi, 4), "]\n\n")

  # --- Bland-Altman data export (CC reconstructs the plot) ------------------
  ba_data <- hosp_rates %>%
    mutate(
      trial          = trial_label,
      mean_rate      = (rate_alg + rate_flowsheet) / 2,
      diff_rate      = rate_alg - rate_flowsheet,
      mean_diff_line = ccc_out$mean_diff,
      loa_lo         = ccc_out$loa_lo,
      loa_hi         = ccc_out$loa_hi
    )

  list(ccc = ccc_out, ba_data = ba_data, hosp_rates = hosp_rates)
}

a6_ccc_sat <- NULL
a6_ccc_sbt <- NULL

if (site_has_flowsheet_sat) {
  a6_ccc_sat <- run_a6_ccc("SAT", "SAT_delivered_primary",
                             "flowsheet_SAT", "SAT_eligible")
} else {
  cat("  SAT CCC/Bland-Altman skipped: no flowsheet SAT data at this site.\n\n")
}

if (site_has_flowsheet_sbt) {
  a6_ccc_sbt <- run_a6_ccc("SBT", "SBT_delivered_2min",
                             "flowsheet_SBT", "SBT_eligible")
} else {
  cat("  SBT CCC/Bland-Altman skipped: no flowsheet SBT data at this site.\n\n")
}

# ---------------------------------------------------------------------------
# 6.5 VISUALIZATIONS: CATERPILLAR + FUNNEL PLOTS
# ---------------------------------------------------------------------------

cat("-- 6.5 Visualizations: caterpillar and funnel plots\n\n")

build_a6_figures <- function(trial_label, blup_df, re_stats) {

  if (single_hospital || all(is.na(blup_df$risk_adj_rate))) {
    cat("  SKIP figures:", trial_label,
        "-- single-hospital site, no hospital-level estimates.\n\n")
    return(NULL)
  }

  blup_plot <- blup_df %>%
    arrange(risk_adj_rate) %>%
    mutate(rank      = row_number(),
           hosp_label = paste0("H", rank))   # de-identified rank label

  grand_mean_rate <- plogis(fixef(
    if (trial_label == "SAT") a6_sat$fit else a6_sbt$fit
  )[["(Intercept)"]])

  # --- Caterpillar plot (main figure) ----------------------------------------
  fig_cat <- ggplot(blup_plot,
                    aes(x = risk_adj_rate,
                        xmin = risk_adj_rate_lo,
                        xmax = risk_adj_rate_hi,
                        y = reorder(hosp_label, risk_adj_rate))) +
    geom_vline(xintercept = grand_mean_rate, linetype = "dashed",
               color = "gray50", linewidth = 0.6) +
    geom_errorbarh(height = 0.3, linewidth = 0.7,
                   color = if (trial_label == "SAT") clr_sat else clr_sbt) +
    geom_point(size = 3,
               color = if (trial_label == "SAT") clr_sat else clr_sbt) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    labs(
      title    = paste0(trial_label,
                        " Delivery: Risk-Adjusted Hospital Rates (BLUPs)"),
      subtitle = paste0("Empirical Bayes shrinkage estimates | ",
                        "Dashed line = grand mean | ",
                        "MOR = ", round(re_stats$MOR, 3),
                        " | VPC = ", round(re_stats$vpc_icc, 3)),
      x        = "Risk-Adjusted Delivery Rate",
      y        = "Hospital (ranked)",
      caption  = paste0(
        "BLUPs shrink small-hospital estimates toward grand mean ",
        "(Rogers et al. 2013).\n",
        "CIs reflect conditional variance of random effects (postVar).\n",
        "Hospital labels de-identified -- H1 = lowest, H",
        nrow(blup_plot), " = highest rate."
      )
    ) +
    theme_abtrise()

  # --- Funnel plot (supplemental) -------------------------------------------
  # x-axis: hospital volume (n eligible days) -- precision proxy
  # y-axis: raw delivery rate
  # Control limits: grand mean ± 1.96*SE and ± 3*SE
  # SE(rate) = sqrt(p*(1-p)/n) where p = grand mean raw rate

  grand_raw <- mean(blup_df$raw_rate, na.rm = TRUE)

  funnel_limits <- blup_df %>%
    summarise(n_range = list(seq(min(n_eligible), max(n_eligible),
                                  length.out = 200))) %>%
    tidyr::unnest(n_range) %>%
    rename(n = n_range) %>%
    mutate(
      se        = sqrt(grand_raw * (1 - grand_raw) / n),
      lim_95_lo = grand_raw - 1.96 * se,
      lim_95_hi = grand_raw + 1.96 * se,
      lim_99_lo = grand_raw - 3.00 * se,
      lim_99_hi = grand_raw + 3.00 * se
    )

  fig_funnel <- ggplot() +
    geom_ribbon(data = funnel_limits,
                aes(x = n, ymin = lim_99_lo, ymax = lim_99_hi),
                fill = "gray85", alpha = 0.6) +
    geom_ribbon(data = funnel_limits,
                aes(x = n, ymin = lim_95_lo, ymax = lim_95_hi),
                fill = "gray70", alpha = 0.6) +
    geom_hline(yintercept = grand_raw, linetype = "dashed",
               color = "gray40", linewidth = 0.6) +
    geom_point(data = blup_df,
               aes(x = n_eligible, y = raw_rate),
               color = if (trial_label == "SAT") clr_sat else clr_sbt,
               size = 3, alpha = 0.85) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    labs(
      title    = paste0(trial_label,
                        " Delivery: Funnel Plot (Raw Rates by Hospital Volume)"),
      subtitle = "Shaded bands: \u00b11.96 SD (dark) and \u00b13 SD (light) control limits",
      x        = "Hospital Volume (eligible days)",
      y        = "Raw Delivery Rate",
      caption  = paste0(
        "Points outside control limits are statistical outliers.\n",
        "Small hospitals naturally show wider scatter -- compare with\n",
        "caterpillar plot (BLUP-shrunk estimates) for full picture.\n",
        "Spiegelhalter (Stat Med 2005)."
      )
    ) +
    theme_abtrise()

  list(caterpillar = fig_cat, funnel = fig_funnel,
       blup_plot_data = blup_plot,
       funnel_limits  = funnel_limits)
}

figs_sat <- build_a6_figures("SAT", a6_sat$blups, a6_sat$re_stats)
figs_sbt <- build_a6_figures("SBT", a6_sbt$blups, a6_sbt$re_stats)

# ---------------------------------------------------------------------------
# 6.6 EXPORTS
# ---------------------------------------------------------------------------

cat("-- 6.6 Exporting Analysis 6 outputs\n\n")

# --- Track 1: Local GLMM outputs -------------------------------------------

# Fixed effect coefficients
export_csv(a6_sat$tidy, "models/a6", "A6_SAT_glmm_coefs.csv")
export_csv(a6_sbt$tidy, "models/a6", "A6_SBT_glmm_coefs.csv")

# RE variance, VPC/ICC, MOR
export_csv(
  bind_rows(a6_sat$re_stats, a6_sbt$re_stats),
  "models/a6", "A6_re_variance_mor.csv"
)

# BLUPs and risk-adjusted rates
export_csv(
  bind_rows(a6_sat$blups, a6_sbt$blups),
  "models/a6", "A6_blups_risk_adj_rates.csv"
)

# ICC(patient) diagnostics
export_csv(
  bind_rows(a6_sat$icc_patient, a6_sbt$icc_patient),
  "models/a6", "A6_icc_patient_diagnostic.csv"
)

# RDS model objects
export_rds(a6_sat$fit, "models/a6", "A6_fit_SAT_glmm.rds")
export_rds(a6_sbt$fit, "models/a6", "A6_fit_SBT_glmm.rds")

# --- Track 2: Hospital aggregate summary (CC pooled model input) -----------
export_csv(hosp_summary_all, "models/a6", "A6_hospital_aggregate_summary.csv")

# --- CCC + Bland-Altman ----------------------------------------------------
if (!is.null(a6_ccc_sat)) {
  export_csv(bind_rows(a6_ccc_sat$ccc),
             "models/a6", "A6_SAT_ccc_results.csv")
  export_csv(a6_ccc_sat$ba_data,
             "figures/a6", "A6_SAT_bland_altman_data.csv")
  export_csv(a6_ccc_sat$hosp_rates,
             "figures/a6", "A6_SAT_hosp_rates_comparison.csv")
}

if (!is.null(a6_ccc_sbt)) {
  export_csv(bind_rows(a6_ccc_sbt$ccc),
             "models/a6", "A6_SBT_ccc_results.csv")
  export_csv(a6_ccc_sbt$ba_data,
             "figures/a6", "A6_SBT_bland_altman_data.csv")
  export_csv(a6_ccc_sbt$hosp_rates,
             "figures/a6", "A6_SBT_hosp_rates_comparison.csv")
}

# --- Figure exports and figure data CSVs for CC reconstruction -------------

if (!is.null(figs_sat)) {
  export_png(figs_sat$caterpillar, "figures/a6",
             "fig_A6_SAT_caterpillar.png", width = 9, height = 7)
  export_png(figs_sat$funnel, "figures/a6",
             "fig_A6_SAT_funnel.png", width = 8, height = 6)
  export_csv(figs_sat$blup_plot_data, "figures/a6",
             "A6_SAT_caterpillar_data.csv")
  export_csv(figs_sat$funnel_limits,  "figures/a6",
             "A6_SAT_funnel_limits_data.csv")
}

if (!is.null(figs_sbt)) {
  export_png(figs_sbt$caterpillar, "figures/a6",
             "fig_A6_SBT_caterpillar.png", width = 9, height = 7)
  export_png(figs_sbt$funnel, "figures/a6",
             "fig_A6_SBT_funnel.png", width = 8, height = 6)
  export_csv(figs_sbt$blup_plot_data, "figures/a6",
             "A6_SBT_caterpillar_data.csv")
  export_csv(figs_sbt$funnel_limits,  "figures/a6",
             "A6_SBT_funnel_limits_data.csv")
}

# Raw delivery rates by hospital -- for CC figure overlay if needed
export_csv(
  bind_rows(
    df_a6_sat %>%
      group_by(hospital_id) %>%
      summarise(trial = "SAT", n_eligible = n(),
                n_delivered = sum(delivered),
                raw_rate = mean(delivered), .groups = "drop"),
    df_a6_sbt %>%
      group_by(hospital_id) %>%
      summarise(trial = "SBT", n_eligible = n(),
                n_delivered = sum(delivered),
                raw_rate = mean(delivered), .groups = "drop")
  ),
  "figures/a6", "A6_raw_rates_by_hospital.csv"
)

# --- Updated waterfall with A6 complete case steps -------------------------
export_csv(waterfall, "diagnostics", "exclusion_waterfall.csv")

cat("\nAnalysis 6 complete.\n\n")


# =============================================================================
# CONSORT / STUDY FLOW DIAGRAMS
# =============================================================================
# Built here because all waterfall steps are now complete (setup, A3-A5, A6).
# Produces two figures -- one per file -- using the consort package.
#
# NOTE ON EXPORT: consort::plot() returns a grid object, not a ggplot.
# Must use png() / dev.off() rather than ggsave() or export_png().
#
# Required package: consort (CRAN). Install with:
#   install.packages("consort")
# =============================================================================

cat("-- Building CONSORT diagrams\n")

if (!requireNamespace("consort", quietly = TRUE)) {
  cat("  NOTE: consort package not installed. Skipping CONSORT diagrams.\n")
  cat("  Install with: install.packages('consort')\n\n")
} else {

  library(consort)

  wf <- waterfall

  # ── Helper: safely pull values from waterfall ────────────────────────────
  pull_n <- function(step_name, file_name) {
    x <- wf %>%
      filter(step == step_name, file == file_name) %>%
      pull(n_remaining)
    if (length(x) == 0) return(NA_integer_)
    as.integer(x[1])
  }

  pull_excl <- function(step_name, file_name) {
    x <- wf %>%
      filter(step == step_name, file == file_name) %>%
      pull(n_excluded)
    if (length(x) == 0) return(NA_integer_)
    as.integer(x[1])
  }

  # ── FILE 1 NUMBERS ────────────────────────────────────────────────────────
  f1_raw        <- pull_n("1_raw_load",                        "File1")
  f1_n_imp_excl <- pull_excl("3_data_coding_error_exclusion",  "File1")
  f1_post_imp   <- pull_n("3_data_coding_error_exclusion",     "File1")
  f1_n_vd_excl  <- pull_excl("3_vent_day_filter",              "File1")
  f1_n_a3_excl  <- pull_excl("5_complete_case_A3",             "File1")
  f1_post_a3    <- pull_n("5_complete_case_A3",                "File1")
  f1_n_a6s_excl <- pull_excl("6_complete_case_A6_SAT",         "File1")
  f1_post_a6s   <- pull_n("6_complete_case_A6_SAT",            "File1")
  f1_n_a6b_excl <- pull_excl("6_complete_case_A6_SBT",         "File1")
  f1_post_a6b   <- pull_n("6_complete_case_A6_SBT",            "File1")

  # ── FILE 2 NUMBERS ────────────────────────────────────────────────────────
  f2_raw_all    <- pull_n("1_raw_load",                        "File2")
  f2_n_zv_excl  <- pull_excl("2_zero_vent_exclusion",          "File2")
  f2_post_zv    <- pull_n("2_zero_vent_exclusion",             "File2")
  f2_n_imp_excl <- pull_excl("3_data_coding_error_exclusion",  "File2")
  f2_post_imp   <- pull_n("3_data_coding_error_exclusion",     "File2")
  f2_n_a4_excl  <- pull_excl("5_complete_case_A4",             "File2")
  f2_post_a4    <- pull_n("5_complete_case_A4",                "File2")
  f2_n_a5_excl  <- pull_excl("5_complete_case_A5",             "File2")
  f2_post_a5    <- pull_n("5_complete_case_A5",                "File2")

  # ── FIGURE 1: FILE 1 (Person-Period / Long) ───────────────────────────────
  cat("  Building File 1 CONSORT diagram...\n")

  tryCatch({

    out_path_f1 <- file.path(
      out_dir, "figures", prefix_file("fig_consort_file1.png")
    )

    png(out_path_f1, width = 10, height = 10, units = "in", res = 300)

    plot(
      add_box(txt = sprintf(
        "%s Person-Days Loaded\n(%s Hospitalizations)\nFile 1 (Person-Period): Raw Data",
        format(f1_raw, big.mark = ","),
        format(n_distinct(df_pp_raw$hospitalization_id), big.mark = ",")
      )) %>%

      add_side_box(txt = sprintf(
        "%s Person-Days Excluded\nImpossible variable values\n(SOFA >24, FiO2 out of range,\nsedation_prior not 0/1)",
        format(f1_n_imp_excl, big.mark = ",")
      ), text_width = 38) %>%

      add_box(txt = sprintf(
        "%s Person-Days\nAfter Data Coding Error Exclusion",
        format(f1_post_imp, big.mark = ",")
      )) %>%

      add_side_box(txt = sprintf(
        "%s Person-Days Excluded\nVentilator day outside 1-28 window",
        format(f1_n_vd_excl, big.mark = ",")
      ), text_width = 38) %>%

      add_box(txt = sprintf(
        "%s Person-Days | %s Hospitalizations\nAnalytic File 1 (vent_day 1-28)",
        format(nrow(df_pp), big.mark = ","),
        format(n_distinct(df_pp$hospitalization_id), big.mark = ",")
      )) %>%

      add_split(txt = c(
        sprintf(
          "Analysis 3\n(Time to Extubation)\n%s Person-Days\n%s Hospitalizations\n(%s excluded: incomplete covariates)",
          format(f1_post_a3,  big.mark = ","),
          format(n_distinct(df_dt$hospitalization_id), big.mark = ","),
          format(f1_n_a3_excl, big.mark = ",")
        ),
        sprintf(
          "Analysis 6\nSAT-Eligible Days\n%s Person-Days\n(%s excluded: incomplete covariates)",
          format(f1_post_a6s,  big.mark = ","),
          format(f1_n_a6s_excl, big.mark = ",")
        ),
        sprintf(
          "Analysis 6\nSBT-Eligible Days\n%s Person-Days\n(%s excluded: incomplete covariates)",
          format(f1_post_a6b,  big.mark = ","),
          format(f1_n_a6b_excl, big.mark = ",")
        )
      ))
    )

    dev.off()
    cat("  File 1 CONSORT diagram exported:", out_path_f1, "\n")

  }, error = function(e) {
    dev.off()
    cat("  WARNING: File 1 CONSORT diagram failed:", conditionMessage(e), "\n")
    cat("  exclusion_waterfall.csv remains available as text alternative.\n")
  })

  # ── FIGURE 2: FILE 2 (Hospitalization / Wide) ─────────────────────────────
  cat("  Building File 2 CONSORT diagram...\n")

  tryCatch({

    out_path_f2 <- file.path(
      out_dir, "figures", prefix_file("fig_consort_file2.png")
    )

    png(out_path_f2, width = 9, height = 9, units = "in", res = 300)

    plot(
      add_box(txt = sprintf(
        "%s Hospitalizations Identified\nFile 2 (Hospitalization-Level): Raw Data",
        format(f2_raw_all, big.mark = ",")
      )) %>%

      add_side_box(txt = sprintf(
        "%s Hospitalizations Excluded\nZero algorithm-captured vent days\n(Fast-extubation survivors,\nmedian VFD = 28)",
        format(f2_n_zv_excl, big.mark = ",")
      ), text_width = 38) %>%

      add_box(txt = sprintf(
        "%s Hospitalizations\nAfter Zero Vent-Day Exclusion",
        format(f2_post_zv, big.mark = ",")
      )) %>%

      add_side_box(txt = sprintf(
        "%s Hospitalizations Excluded\nImpossible variable values\n(age <18, CCI >37, VFD_28 >28,\nSOFA/FiO2/PEEP out of range)",
        format(f2_n_imp_excl, big.mark = ",")
      ), text_width = 38) %>%

      add_box(txt = sprintf(
        "%s Hospitalizations\nAnalytic Cohort (File 2)",
        format(f2_post_imp, big.mark = ",")
      )) %>%

      add_split(txt = c(
        sprintf(
          "Analysis 4\n(VFD-28 Two-Part Model)\n%s Hospitalizations\n(%s excluded: incomplete covariates)",
          format(f2_post_a4,  big.mark = ","),
          format(f2_n_a4_excl, big.mark = ",")
        ),
        sprintf(
          "Analysis 5\n(ICU LOS / Mortality)\n%s Hospitalizations\n(%s excluded: incomplete covariates\nor ICU LOS <1)",
          format(f2_post_a5,  big.mark = ","),
          format(f2_n_a5_excl, big.mark = ",")
        )
      ))
    )

    dev.off()
    cat("  File 2 CONSORT diagram exported:", out_path_f2, "\n\n")

  }, error = function(e) {
    dev.off()
    cat("  WARNING: File 2 CONSORT diagram failed:", conditionMessage(e), "\n")
    cat("  exclusion_waterfall.csv remains available as text alternative.\n\n")
  })

} # end consort block

# =============================================================================
# SESSION INFO
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "diagnostics", prefix_file("session_info_a6.txt")))
cat("Session info saved.\n")
cat("\n=== Analysis 6 script complete ===\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
