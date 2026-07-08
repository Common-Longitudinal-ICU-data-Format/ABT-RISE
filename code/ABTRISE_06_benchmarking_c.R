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
#   outputs/A6_benchmark_outcomes/models/   A6_SAT_glmm_coefs.csv, A6_SBT_glmm_coefs.csv,
#                                          A6_re_variance_mor.csv, A6_blups_risk_adj_rates.csv,
#                                          A6_icc_patient_diagnostic.csv,
#                                          A6_hospital_aggregate_summary.csv,
#                                          A6_SAT_ccc_results.csv, A6_SBT_ccc_results.csv [if flowsheet]
#   outputs/A6_benchmark_outcomes/figures/  fig_A6_SAT_caterpillar.png, fig_A6_SBT_caterpillar.png,
#                                          fig_A6_SAT_funnel.png, fig_A6_SBT_funnel.png,
#                                          fig_A6_SAT_ranked_hospitals.png, fig_A6_SBT_ranked_hospitals.png,
#                                          A6_*_caterpillar_data.csv, A6_*_funnel_limits_data.csv,
#                                          A6_*_ranked_hospitals_data.csv,
#                                          A6_*_bland_altman_data.csv [if flowsheet]
#   outputs/diagnostics/                   fig_consort_file1.png, fig_consort_file2.png,
#                                          session_info_a6.txt, exclusion_waterfall.csv
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
#                      "sedation_prior", "NEE_prior")

# --- A6 SAT dataset: eligible SAT days, complete case ----------------------
df_a6_sat <- df_pp %>%
  filter(SAT_eligible == 1L) %>%
  filter(
    !is.na(SAT_delivered_primary),
    !is.na(age), !is.na(sex), !is.na(CCI),
    !is.na(SOFA_prior), !is.na(FiO2_prior), !is.na(PEEP_prior),
    !is.na(sedation_prior),
  
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
     # hospital_type      = first(as.character(hospital_type)),
      #location_type      = first(as.character(location_type)),
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
    n_hospitals    = as.integer(n_hosp_fs),
    CCC            = as.numeric(ccc_result$rho.c$est[1]),    # as.numeric([1]): force scalar;
    CCC_lo         = as.numeric(ccc_result$rho.c$lower[1]),  # epiR version differences can return
    CCC_hi         = as.numeric(ccc_result$rho.c$upper[1]),  # named numeric, 1-row df, or list --
    pearson_r      = as.numeric(ccc_result$r[1]),            # any of these stored in a tibble col
    cb_coef        = as.numeric(ccc_result$C.b[1]),          # becomes a list col write_csv rejects
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
        nrow(blup_plot), " = highest rate.\n"
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
        "Spiegelhalter (Stat Med 2005).\n"
      )
    ) +
    theme_abtrise()

  # --- Hospital ranking figure (vertical caterpillar) -------------------------
  # Style: ranked-hospital dot plot, colored by hospital type.
  # Mirrors: Figure 3, Valk et al., Critical Care Medicine (proning study).
  # X-axis: hospitals ranked by raw delivery rate (integer rank labels hidden;
  #         de-identified to H1...Hn). Y-axis: raw delivery rate with Wilson CI.
  # Color: hospital_type (academic / community / other) from df_pp.
  # Gated on multi-hospital sites and n_delivered > 0 (CI undefined at n=0).

  # Wilson CI helper (base R -- no new packages)
  wilson_ci <- function(x, n, conf = 0.95) {
    z  <- qnorm(1 - (1 - conf) / 2)
    p  <- x / n
    cn <- (p + z^2 / (2*n)) / (1 + z^2/n)
    mg <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n)
    list(lo = pmax(0, cn - mg), hi = pmin(1, cn + mg))
  }

  # Join hospital_type from df_pp (available in global scope from setup)
  hosp_type_lkp <- df_pp %>%
    group_by(hospital_id) %>%
    summarise(hospital_type = tolower(first(as.character(hospital_type))),
              .groups = "drop")

  ranking_data <- blup_df %>%
    left_join(hosp_type_lkp, by = "hospital_id") %>%
    mutate(
      hospital_type = replace_na(hospital_type, "unknown"),
      rank          = rank(raw_rate, ties.method = "first"),
      hosp_label    = paste0("H", rank)
    ) %>%
    filter(!is.na(n_delivered), n_eligible > 0) %>%
    mutate(
      ci_lo = wilson_ci(n_delivered, n_eligible)$lo,
      ci_hi = wilson_ci(n_delivered, n_eligible)$hi
    )

  # Hospital type color scale: reuse JAMA palette subsets without touching
  # clr_sat / clr_sbt (trial colors). Academic=red, Community=blue, Other=gray.
  hosp_type_vals <- c(
    "academic"  = JAMA_COLORS[4],  # muted red
    "community" = JAMA_COLORS[3],  # JAMA blue
    "unknown"   = JAMA_COLORS[7],  # warm gray
    "other"     = JAMA_COLORS[7]
  )

  fig_ranking <- ggplot(
      ranking_data,
      aes(x = rank, y = raw_rate, ymin = ci_lo, ymax = ci_hi,
          color = hospital_type)
    ) +
    geom_hline(yintercept = mean(ranking_data$raw_rate, na.rm = TRUE),
               linetype = "dashed", color = "gray45", linewidth = 0.55) +
    geom_errorbar(width = 0.3, linewidth = 0.7, alpha = 0.6) +
    geom_point(size = 3, alpha = 0.9) +
    scale_x_continuous(
      breaks  = ranking_data$rank,
      labels  = ranking_data$hosp_label,
      expand  = expansion(add = 0.8)
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.04))
    ) +
    scale_color_manual(
      values = hosp_type_vals,
      name   = "Hospital Type",
      na.value = JAMA_COLORS[7]
    ) +
    labs(
      title    = paste0(trial_label,
                        " Delivery: Hospitals Ranked by Raw Delivery Rate"),
      subtitle = paste0(
        "Error bars = 95% Wilson CI  |  Dashed line = consortium mean  |  ",
        "n = ", nrow(ranking_data), " hospitals"
      ),
      x        = paste0("Hospitals\nRanked by ", trial_label,
                        " Raw Delivery Rate"),
      y        = paste0("% Eligible Days With ", trial_label, " Delivered"),
      caption  = paste0(
        "Hospital labels de-identified (H1 = lowest, H",
        nrow(ranking_data), " = highest raw rate).\n",
        "Wilson 95% CI from raw eligible-day counts.\n",
        "Hospital type from site data -- may be unlabeled if not captured.\n",
        "Reference: Valk et al. (proning study figure style)."
      )
    ) +
    theme_abtrise() +
    theme(
      legend.position  = "top",
      axis.text.x      = element_text(size = 7, angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  list(caterpillar    = fig_cat,
       funnel         = fig_funnel,
       ranking        = fig_ranking,
       blup_plot_data = blup_plot,
       funnel_limits  = funnel_limits,
       ranking_data   = ranking_data)
}

figs_sat <- build_a6_figures("SAT", a6_sat$blups, a6_sat$re_stats)
figs_sbt <- build_a6_figures("SBT", a6_sbt$blups, a6_sbt$re_stats)

# ---------------------------------------------------------------------------
# 6.6 EXPORTS
# ---------------------------------------------------------------------------

cat("-- 6.6 Exporting Analysis 6 outputs\n\n")

# --- Track 1: Local GLMM outputs -------------------------------------------

# Fixed effect coefficients
export_csv(a6_sat$tidy, "A6_benchmark_outcomes/models", "A6_SAT_glmm_coefs.csv")
export_csv(a6_sbt$tidy, "A6_benchmark_outcomes/models", "A6_SBT_glmm_coefs.csv")

# RE variance, VPC/ICC, MOR
export_csv(
  bind_rows(a6_sat$re_stats, a6_sbt$re_stats),
  "A6_benchmark_outcomes/models", "A6_re_variance_mor.csv"
)

# BLUPs and risk-adjusted rates
export_csv(
  bind_rows(a6_sat$blups, a6_sbt$blups),
  "A6_benchmark_outcomes/models", "A6_blups_risk_adj_rates.csv"
)

# ICC(patient) diagnostics
export_csv(
  bind_rows(a6_sat$icc_patient, a6_sbt$icc_patient),
  "A6_benchmark_outcomes/models", "A6_icc_patient_diagnostic.csv"
)

# --- Track 2: Hospital aggregate summary (CC pooled model input) -----------
export_csv(hosp_summary_all, "A6_benchmark_outcomes/models", "A6_hospital_aggregate_summary.csv")

# SA: Hospital aggregate summary split by age group (<65 vs >=65)
hosp_summary_age_split <- bind_rows(
  build_hosp_summary("SAT", df_a6_sat %>% filter(age <  65)) %>% mutate(age_group = "age_u65"),
  build_hosp_summary("SAT", df_a6_sat %>% filter(age >= 65)) %>% mutate(age_group = "age_ge65"),
  build_hosp_summary("SBT", df_a6_sbt %>% filter(age <  65)) %>% mutate(age_group = "age_u65"),
  build_hosp_summary("SBT", df_a6_sbt %>% filter(age >= 65)) %>% mutate(age_group = "age_ge65")
)
export_csv(hosp_summary_age_split, "A6_benchmark_outcomes/models", "SA_age_split_A6_hospital_aggregate_summary.csv")

# --- CCC + Bland-Altman ----------------------------------------------------
if (!is.null(a6_ccc_sat)) {
  export_csv(bind_rows(a6_ccc_sat$ccc),
             "A6_benchmark_outcomes/models", "A6_SAT_ccc_results.csv")
  export_csv(a6_ccc_sat$ba_data,
             "A6_benchmark_outcomes/figures", "A6_SAT_bland_altman_data.csv")
  export_csv(a6_ccc_sat$hosp_rates,
             "A6_benchmark_outcomes/figures", "A6_SAT_hosp_rates_comparison.csv")
}

if (!is.null(a6_ccc_sbt)) {
  export_csv(bind_rows(a6_ccc_sbt$ccc),
             "A6_benchmark_outcomes/models", "A6_SBT_ccc_results.csv")
  export_csv(a6_ccc_sbt$ba_data,
             "A6_benchmark_outcomes/figures", "A6_SBT_bland_altman_data.csv")
  export_csv(a6_ccc_sbt$hosp_rates,
             "A6_benchmark_outcomes/figures", "A6_SBT_hosp_rates_comparison.csv")
}

# --- Figure exports and figure data CSVs for CC reconstruction -------------

if (!is.null(figs_sat)) {
  export_png(figs_sat$caterpillar, "A6_benchmark_outcomes/figures",
             "fig_A6_SAT_caterpillar.png", width = 9, height = 7)
  export_png(figs_sat$funnel, "A6_benchmark_outcomes/figures",
             "fig_A6_SAT_funnel.png", width = 8, height = 6)
  export_png(figs_sat$ranking, "A6_benchmark_outcomes/figures",
             "fig_A6_SAT_ranked_hospitals.png", width = 10, height = 6)
  export_csv(figs_sat$blup_plot_data, "A6_benchmark_outcomes/figures",
             "A6_SAT_caterpillar_data.csv")
  export_csv(figs_sat$funnel_limits,  "A6_benchmark_outcomes/figures",
             "A6_SAT_funnel_limits_data.csv")
  export_csv(figs_sat$ranking_data,   "A6_benchmark_outcomes/figures",
             "A6_SAT_ranked_hospitals_data.csv")
}

if (!is.null(figs_sbt)) {
  export_png(figs_sbt$caterpillar, "A6_benchmark_outcomes/figures",
             "fig_A6_SBT_caterpillar.png", width = 9, height = 7)
  export_png(figs_sbt$funnel, "A6_benchmark_outcomes/figures",
             "fig_A6_SBT_funnel.png", width = 8, height = 6)
  export_png(figs_sbt$ranking, "A6_benchmark_outcomes/figures",
             "fig_A6_SBT_ranked_hospitals.png", width = 10, height = 6)
  export_csv(figs_sbt$blup_plot_data, "A6_benchmark_outcomes/figures",
             "A6_SBT_caterpillar_data.csv")
  export_csv(figs_sbt$funnel_limits,  "A6_benchmark_outcomes/figures",
             "A6_SBT_funnel_limits_data.csv")
  export_csv(figs_sbt$ranking_data,   "A6_benchmark_outcomes/figures",
             "A6_SBT_ranked_hospitals_data.csv")
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
  "A6_benchmark_outcomes/figures", "A6_raw_rates_by_hospital.csv"
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

# =============================================================================
# CONSORT / STUDY FLOW DIAGRAMS (ggplot2 implementation)
# =============================================================================
# Built here because all waterfall steps are now complete.
# Style: main flow boxes (white/light, dark border) + exclusion side bars
# (JAMA blue filled, white text), placed left of the main flow.
# Mirrors the visual convention of: Rotta et al., BMC Health Services Research
# 2022 (DOI:10.1186/s12913-022-07467-8), Figure 1.
#
# Replaces the consort-package approach (removed) with a self-contained
# ggplot2 function. No additional packages required beyond those already loaded.
#
# draw_consort_gg():
#   steps: named list of lists, each with:
#     $main_txt:  character -- main box text
#     $excl_txt:  character or NULL -- exclusion left-bar text
#     $split_txts: character vector or NULL -- if provided, the final main box
#                  splits into N equal-width child boxes at the bottom
# =============================================================================

cat("-- Building CONSORT diagrams (ggplot2, blue exclusion bars)\n")

draw_consort_gg <- function(steps, title_txt = NULL) {

  # --- Layout constants -------------------------------------------------------
  # Main flow: centered at x = 0.62; exclusion bars: centered at x = 0.18
  main_cx  <- 0.62;  main_bw  <- 0.40;  main_bh  <- 0.085
  excl_cx  <- 0.18;  excl_bw  <- 0.30;  excl_bh  <- 0.070
  arrow_col <- "gray40"

  n_main_steps <- sum(sapply(steps, function(s) is.null(s$split_txts)))
  n_total      <- length(steps)
  has_split    <- !is.null(steps[[n_total]]$split_txts)

  # Y positions for main boxes (top to bottom), leaving room for split row
  split_y_top <- 0.09
  y_tops <- seq(0.96,
                ifelse(has_split, split_y_top + main_bh + 0.05, 0.09),
                length.out = ifelse(has_split, n_total - 1, n_total))

  # --- Accumulate geom data --------------------------------------------------
  df_main   <- tibble()
  df_excl   <- tibble()
  df_segs   <- tibble()
  df_split  <- tibble()

  for (i in seq_along(y_tops)) {
    step <- steps[[i]]
    cy   <- y_tops[i] - main_bh / 2

    df_main <- bind_rows(df_main,
      tibble(xmin = main_cx - main_bw/2, xmax = main_cx + main_bw/2,
             ymin = cy - main_bh/2,      ymax = cy + main_bh/2,
             cx = main_cx, cy = cy, txt = step$main_txt))

    # Downward connector to next step (or to split bar)
    next_top_y <- if (i < length(y_tops)) {
      y_tops[i+1] - main_bh / 2 + main_bh / 2   # top of next main box
    } else if (has_split) {
      split_y_top + main_bh   # top of split row
    } else NULL

    if (!is.null(next_top_y)) {
      df_segs <- bind_rows(df_segs,
        tibble(x = main_cx, xend = main_cx,
               y = cy - main_bh/2,       # bottom of this box
               yend = next_top_y,
               seg_type = "main_flow"))
    }

    # Exclusion side bar + horizontal connector
    if (!is.null(step$excl_txt)) {
      excl_cy <- cy
      df_excl <- bind_rows(df_excl,
        tibble(xmin = excl_cx - excl_bw/2, xmax = excl_cx + excl_bw/2,
               ymin = excl_cy - excl_bh/2, ymax = excl_cy + excl_bh/2,
               cx = excl_cx, cy = excl_cy, txt = step$excl_txt))
      # Horizontal connector: left edge of main box → right edge of excl bar
      df_segs <- bind_rows(df_segs,
        tibble(x = main_cx - main_bw/2, xend = excl_cx + excl_bw/2,
               y = excl_cy,             yend = excl_cy,
               seg_type = "excl_horiz"))
    }
  }

  # --- Split row --------------------------------------------------------------
  if (has_split) {
    split_txts <- steps[[n_total]]$split_txts
    n_split    <- length(split_txts)
    split_bw   <- min(0.30, 0.80 / n_split)
    split_cy   <- split_y_top + main_bh / 2
    split_xs   <- seq(0.22, 0.96, length.out = n_split)

    for (j in seq_len(n_split)) {
      df_split <- bind_rows(df_split,
        tibble(xmin = split_xs[j] - split_bw/2,
               xmax = split_xs[j] + split_bw/2,
               ymin = split_cy - main_bh/2,
               ymax = split_cy + main_bh/2,
               cx   = split_xs[j], cy = split_cy,
               txt  = split_txts[j]))
    }

    # Horizontal span bar connecting all split boxes
    span_y <- split_y_top + main_bh + 0.02
    df_segs <- bind_rows(df_segs,
      tibble(x = split_xs[1], xend = split_xs[n_split],
             y = span_y, yend = span_y, seg_type = "split_span"))
    # Vertical drops from span bar to each split box top
    for (j in seq_len(n_split)) {
      df_segs <- bind_rows(df_segs,
        tibble(x = split_xs[j], xend = split_xs[j],
               y = span_y, yend = split_cy + main_bh/2,
               seg_type = "split_drop"))
    }
  }

  # --- Assemble ggplot --------------------------------------------------------
  p <- ggplot() +
    # Connectors (drawn first, behind boxes)
    geom_segment(data = df_segs,
                 aes(x=x, xend=xend, y=y, yend=yend),
                 color = arrow_col, linewidth = 0.55) +
    # Main flow boxes: white fill, gray border
    geom_rect(data = df_main,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = "white", color = "gray35", linewidth = 0.45) +
    geom_text(data = df_main,
              aes(x=cx, y=cy, label=txt),
              size = 2.65, color = "gray15", hjust = 0.5, vjust = 0.5,
              lineheight = 1.05) +
    # Exclusion bars: JAMA blue fill, white text (BMC Health Serv Res style)
    {if (nrow(df_excl) > 0)
        list(
          geom_rect(data = df_excl,
                    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                    fill = JAMA_COLORS[3], color = JAMA_COLORS[3],
                    linewidth = 0.4),
          geom_text(data = df_excl,
                    aes(x=cx, y=cy, label=txt),
                    size = 2.3, color = "white", hjust = 0.5, vjust = 0.5,
                    lineheight = 0.95)
        )
    } +
    # Split boxes (bottom row): same style as main boxes
    {if (nrow(df_split) > 0)
        list(
          geom_rect(data = df_split,
                    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                    fill = "white", color = "gray35", linewidth = 0.45),
          geom_text(data = df_split,
                    aes(x=cx, y=cy, label=txt),
                    size = 2.3, color = "gray15", hjust = 0.5, vjust = 0.5,
                    lineheight = 0.95)
        )
    } +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()

  if (!is.null(title_txt))
    p <- p + labs(title = title_txt) +
             theme(plot.title = element_text(size = 9, hjust = 0.5,
                                             margin = margin(b = 4)))
  p
}

# ── Pull waterfall values (same helper pattern as before) -------------------
pull_n    <- function(step_name, file_name) {
  x <- waterfall %>% filter(step == step_name, file == file_name) %>%
    pull(n_remaining)
  if (length(x) == 0) return(NA_integer_)
  as.integer(x[1])
}
pull_excl <- function(step_name, file_name) {
  x <- waterfall %>% filter(step == step_name, file == file_name) %>%
    pull(n_excluded)
  if (length(x) == 0) return(NA_integer_)
  as.integer(x[1])
}

f1_raw        <- pull_n   ("1_raw_load",                       "File1")
f1_n_imp_excl <- pull_excl("3_data_coding_error_exclusion",    "File1")
f1_post_imp   <- pull_n   ("3_data_coding_error_exclusion",    "File1")
f1_n_vd_excl  <- pull_excl("3_vent_day_filter",               "File1")
f1_n_a3_excl  <- pull_excl("5_complete_case_A3",              "File1")
f1_post_a3    <- pull_n   ("5_complete_case_A3",              "File1")
f1_n_a6s_excl <- pull_excl("6_complete_case_A6_SAT",          "File1")
f1_post_a6s   <- pull_n   ("6_complete_case_A6_SAT",          "File1")
f1_n_a6b_excl <- pull_excl("6_complete_case_A6_SBT",          "File1")
f1_post_a6b   <- pull_n   ("6_complete_case_A6_SBT",          "File1")

f2_raw_all    <- pull_n   ("1_raw_load",                       "File2")
f2_n_zv_excl  <- pull_excl("2_zero_vent_exclusion",           "File2")
f2_post_zv    <- pull_n   ("2_zero_vent_exclusion",           "File2")
f2_n_imp_excl <- pull_excl("3_data_coding_error_exclusion",   "File2")
f2_post_imp   <- pull_n   ("3_data_coding_error_exclusion",   "File2")
f2_n_a4_excl  <- pull_excl("5_complete_case_A4",             "File2")
f2_post_a4    <- pull_n   ("5_complete_case_A4",             "File2")
f2_n_a5_excl  <- pull_excl("5_complete_case_A5",             "File2")
f2_post_a5    <- pull_n   ("5_complete_case_A5",             "File2")

# ── FILE 1 CONSORT DIAGRAM ──────────────────────────────────────────────────
cat("  Building File 1 CONSORT diagram (ggplot2)...\n")

steps_f1 <- list(
  list(
    main_txt  = sprintf(
      "%s Person-Days Loaded\n(%s Hospitalizations)\nFile 1 (Person-Period) — Raw Data",
      format(f1_raw, big.mark = ","),
      format(n_distinct(df_pp_raw$hospitalization_id), big.mark = ",")
    ),
    excl_txt  = NULL,
    split_txts = NULL
  ),
  list(
    main_txt  = sprintf(
      "%s Person-Days\nAfter Data Coding Error Exclusion",
      format(f1_post_imp, big.mark = ",")
    ),
    excl_txt  = sprintf(
      "%s Person-Days Excluded\nImpossible values:\nSOFA >24, FiO2 out of range\nor sedation_prior not 0/1",
      format(f1_n_imp_excl, big.mark = ",")
    ),
    split_txts = NULL
  ),
  list(
    main_txt  = sprintf(
      "%s Person-Days | %s Hospitalizations\nAnalytic File 1 (vent_day 1–28)",
      format(nrow(df_pp), big.mark = ","),
      format(n_distinct(df_pp$hospitalization_id), big.mark = ",")
    ),
    excl_txt  = sprintf(
      "%s Person-Days Excluded\nVent day outside 1–28 window",
      format(f1_n_vd_excl, big.mark = ",")
    ),
    split_txts = NULL
  ),
  list(
    main_txt   = NULL,
    excl_txt   = NULL,
    split_txts = c(
      sprintf(
        "Analysis 3\n(Time to Extubation)\n%s Person-Days\n%s Hospitalizations\n(%s excluded: incomplete\ncovariates)",
        format(f1_post_a3, big.mark = ","),
        format(n_distinct(df_dt$hospitalization_id), big.mark = ","),
        format(f1_n_a3_excl, big.mark = ",")
      ),
      sprintf(
        "Analysis 6\nSAT-Eligible Days\n%s Person-Days\n(%s excluded:\nincomplete covariates)",
        format(f1_post_a6s, big.mark = ","),
        format(f1_n_a6s_excl, big.mark = ",")
      ),
      sprintf(
        "Analysis 6\nSBT-Eligible Days\n%s Person-Days\n(%s excluded:\nincomplete covariates)",
        format(f1_post_a6b, big.mark = ","),
        format(f1_n_a6b_excl, big.mark = ",")
      )
    )
  )
)

fig_consort_f1 <- tryCatch(
  draw_consort_gg(steps_f1, "Study Flow — File 1 (Person-Period)"),
  error = function(e) {
    cat("  WARNING: File 1 CONSORT diagram failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(fig_consort_f1)) {
  export_png(fig_consort_f1, "diagnostics",
             "fig_consort_file1.png", width = 10, height = 10)
  cat("  File 1 CONSORT diagram exported.\n")
} else {
  cat("  exclusion_waterfall.csv remains available as text alternative.\n")
}

# ── FILE 2 CONSORT DIAGRAM ──────────────────────────────────────────────────
cat("  Building File 2 CONSORT diagram (ggplot2)...\n")

steps_f2 <- list(
  list(
    main_txt  = sprintf(
      "%s Hospitalizations Identified\nFile 2 (Hospitalization-Level) — Raw Data",
      format(f2_raw_all, big.mark = ",")
    ),
    excl_txt  = NULL,
    split_txts = NULL
  ),
  list(
    main_txt  = sprintf(
      "%s Hospitalizations\nAfter Zero Vent-Day Exclusion",
      format(f2_post_zv, big.mark = ",")
    ),
    excl_txt  = sprintf(
      "%s Hospitalizations Excluded\nZero algorithm-captured vent days\n(Fast-extubation survivors;\nmedian VFD = 28)",
      format(f2_n_zv_excl, big.mark = ",")
    ),
    split_txts = NULL
  ),
  list(
    main_txt  = sprintf(
      "%s Hospitalizations\nAnalytic Cohort (File 2)",
      format(f2_post_imp, big.mark = ",")
    ),
    excl_txt  = sprintf(
      "%s Hospitalizations Excluded\nImpossible values:\nage <18, CCI >37, VFD_28 >28\nor SOFA/FiO2/PEEP out of range",
      format(f2_n_imp_excl, big.mark = ",")
    ),
    split_txts = NULL
  ),
  list(
    main_txt   = NULL,
    excl_txt   = NULL,
    split_txts = c(
      sprintf(
        "Analysis 4\n(VFD-28 Two-Part Model)\n%s Hospitalizations\n(%s excluded:\nincomplete covariates)",
        format(f2_post_a4, big.mark = ","),
        format(f2_n_a4_excl, big.mark = ",")
      ),
      sprintf(
        "Analysis 5\n(ICU LOS / Mortality)\n%s Hospitalizations\n(%s excluded:\nincomplete covariates\nor ICU LOS < 1)",
        format(f2_post_a5, big.mark = ","),
        format(f2_n_a5_excl, big.mark = ",")
      )
    )
  )
)

fig_consort_f2 <- tryCatch(
  draw_consort_gg(steps_f2, "Study Flow — File 2 (Hospitalization-Level)"),
  error = function(e) {
    cat("  WARNING: File 2 CONSORT diagram failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(fig_consort_f2)) {
  export_png(fig_consort_f2, "diagnostics",
             "fig_consort_file2.png", width = 9, height = 9)
  cat("  File 2 CONSORT diagram exported.\n\n")
} else {
  cat("  exclusion_waterfall.csv remains available as text alternative.\n\n")
}

# =============================================================================
# SESSION INFO
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "diagnostics", prefix_file("session_info_a6.txt")))
cat("Session info saved.\n")
cat("\n=== Analysis 6 script complete ===\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
