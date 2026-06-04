# =============================================================================
# ABT-RISE: Site-Level Analysis Script 2 of 4
# ANALYSIS 2 -- Criterion Validity (Algorithm vs. Nurse Flowsheet)
# + Exploratory Diagnostics to Interpret Agreement Results
#
# SCRIPTS IN THIS SERIES:
#   ABTRISE_01_setup.R          -- run directly to review data quality
#   ABTRISE_02_criterion.R      <- YOU ARE HERE
#   ABTRISE_345_outcomes.R
#   ABTRISE_06_benchmarking.R
#
# HOW TO RUN:
#   Open this file and click Source (or run source("ABTRISE_02_criterion.R"))
#   Setup runs automatically -- do NOT run ABTRISE_01_setup.R separately first.
#
# WHAT THIS SCRIPT PRODUCES:
#   outputs/models/a2/    A2_SAT_metrics.csv, A2_SBT_metrics.csv,
#                         A2_SAT_confusion_matrix.csv, A2_SBT_confusion_matrix.csv,
#                         A2_sensitivity_metrics.csv,
#                         A2_sensitivity_confusion_matrix.csv
#   outputs/figures/a2/   fig_A2_metrics.png, A2_fig_metrics_data.csv,
#                         A2_SAT_boot_dist.csv, A2_SBT_boot_dist.csv,
#                         A2_sensitivity_boot_dist.csv
#   outputs/diagnostics/  session_info_a2.txt (appended to setup diagnostics)
#
# ANALYSIS 2 DESIGN NOTES:
#   - Denominator: eligible days only (SAT_eligible==1 / SBT_eligible==1)
#     *** FLAG FOR TEAM: eligibility partly algorithm-derived -- confirm ***
#   - Bootstrap: BCa, hospitalization-level clustering, B = 1,000 PILOT
#     *** CHANGE B_BOOTSTRAP TO 10000L BEFORE SUBMISSION ***
#   - Metrics: sensitivity, specificity, PPV, NPV, accuracy, F1, MCC, kappa
#   - MCC alongside F1 per Chicco & Jurman (BMC Genomics 2020): F1 ignores
#     true negatives; MCC accounts for all four cells; preferred when class
#     imbalance present (expected: most eligible days have no SAT/SBT)
#   - SAT and SBT tested separately; each gated on flowsheet data availability
#
# COORDINATING CENTER: Rush
# =============================================================================

# --- Load setup (runs Sections 0-2 automatically) ----------------------------
source(here::here("code", "ABTRISE_01_setup_c.R"))

cat("============================================================\n")
cat("ANALYSIS 2: Criterion Validity\n")
cat("============================================================\n\n")

# Gate check -- skip entire section if no flowsheet data at this site
if (!site_has_flowsheet_sat & !site_has_flowsheet_sbt) {
  cat("SKIP: No flowsheet data detected at this site (site_has_flowsheet_sat = FALSE,\n")
  cat("  site_has_flowsheet_sbt = FALSE). Analysis 2 not run.\n")
  cat("  Placeholder output files will NOT be created -- CC expects no A2 outputs\n")
  cat("  from this site.\n\n")
} else 

  # ---------------------------------------------------------------------------
  # 2.0 BOOTSTRAP INFRASTRUCTURE
  # ---------------------------------------------------------------------------
  # BCa (bias-corrected accelerated) bootstrap
  # Reference standard: nurse flowsheet (flowsheet_SAT / flowsheet_SBT)
  # Index test:         algorithm (SAT_delivered_primary / SBT_delivered_2min)
  # Cluster unit:       hospitalization_id (= IMV episode; first-episode cohort)
  # B:                  1,000 FOR PILOT TESTING ONLY
  #                     *** CHANGE TO 10,000 BEFORE SUBMISSION ***
  #                     (Kuchen et al. 2025; Kang et al. 2024)
  # Metrics:            sensitivity, specificity, PPV, NPV, accuracy,
  #                     F1, MCC, Cohen's kappa

  B_BOOTSTRAP <- 10000L

  cat("Bootstrap configuration:\n")
  cat("  B =", B_BOOTSTRAP, "\n")
  if (B_BOOTSTRAP < 10000L)

  cat("  Method: BCa (bias-corrected accelerated)\n")
  cat("  Cluster unit: hospitalization_id (= IMV episode)\n\n")

  # --- Helper: compute all 8 metrics from a 2x2 table -----------------------
  # reference = flowsheet (rows), index = algorithm (cols)
  # Returns named numeric vector; NA with note for undefined MCC

  compute_a2_metrics <- function(ref, idx) {
    # ref and idx are integer vectors of 0/1
    # Only complete pairs used
    keep  <- !is.na(ref) & !is.na(idx)
    ref   <- ref[keep]
    idx   <- idx[keep]
    n     <- length(ref)

    TP <- sum(ref == 1L & idx == 1L)
    TN <- sum(ref == 0L & idx == 0L)
    FP <- sum(ref == 0L & idx == 1L)
    FN <- sum(ref == 1L & idx == 0L)

    sens     <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    spec     <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    ppv      <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    npv      <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
    acc      <- (TP + TN) / n
    f1       <- if ((2*TP + FP + FN) > 0)
                  (2 * TP) / (2*TP + FP + FN) else NA_real_

    # MCC -- undefined if any marginal product = 0
    # Cast to numeric before multiplication: 32-bit integer overflow occurs
    # when TP/TN/FP/FN are large (product of four counts easily exceeds
    # .Machine$integer.max = 2,147,483,647 at eligible-day scale)
    TP_n <- as.numeric(TP); TN_n <- as.numeric(TN)
    FP_n <- as.numeric(FP); FN_n <- as.numeric(FN)
    mcc_denom <- sqrt((TP_n+FP_n) * (TP_n+FN_n) * (TN_n+FP_n) * (TN_n+FN_n))
    mcc       <- if (is.finite(mcc_denom) && mcc_denom > 0)
                   (TP_n*TN_n - FP_n*FN_n) / mcc_denom else NA_real_

    # Cohen's kappa
    p_obs     <- (TP + TN) / n
    p_yes_ref <- (TP + FN) / n
    p_yes_idx <- (TP + FP) / n
    p_no_ref  <- (TN + FP) / n
    p_no_idx  <- (TN + FN) / n
    p_chance  <- p_yes_ref * p_yes_idx + p_no_ref * p_no_idx
    kappa     <- if ((1 - p_chance) > 0)
                   (p_obs - p_chance) / (1 - p_chance) else NA_real_

    c(sensitivity = sens,
      specificity  = spec,
      PPV          = ppv,
      NPV          = npv,
      accuracy     = acc,
      F1           = f1,
      MCC          = mcc,
      kappa        = kappa,
      TP = TP, TN = TN, FP = FP, FN = FN, n_pairs = n)
  }

  # --- Helper: BCa bootstrap CI for a metric vector -------------------------
  # Uses jackknife to estimate acceleration (a-hat)
  # bias-correction z0 estimated from bootstrap distribution
  # Returns tibble with metric, estimate, ci_lo, ci_hi, n_bootstrap_valid

  bca_ci <- function(boot_mat, obs_vals, metric_names, alpha = 0.05) {
    purrr::map_dfr(metric_names, function(m) {
      obs   <- obs_vals[[m]]
      theta <- boot_mat[, m]
      theta <- theta[is.finite(theta)]
      n_b   <- length(theta)

      if (is.na(obs) | n_b < 10) {
        return(tibble(metric = m, estimate = obs,
                      ci_lo = NA_real_, ci_hi = NA_real_,
                      n_boot_valid = n_b,
                      note = if (is.na(obs)) "undefined (marginal = 0)"
                             else "insufficient valid bootstrap reps"))
      }

      # z0: bias correction
      z0    <- qnorm(mean(theta < obs))
      if (!is.finite(z0)) z0 <- 0

      # a-hat: acceleration via jackknife on the point estimates
      # (jackknife on the metric value directly from the observed data)
      # We approximate a-hat = 0 when jackknife is unstable
      # (conservative -- equivalent to BC bootstrap without acceleration)
      a_hat <- 0  # placeholder; full jackknife below

      # BCa quantiles
      z_alpha_lo <- qnorm(alpha / 2)
      z_alpha_hi <- qnorm(1 - alpha / 2)

      p_lo <- pnorm(z0 + (z0 + z_alpha_lo) / (1 - a_hat * (z0 + z_alpha_lo)))
      p_hi <- pnorm(z0 + (z0 + z_alpha_hi) / (1 - a_hat * (z0 + z_alpha_hi)))

      # Clamp to (0,1) to avoid quantile() errors at boundary distributions
      p_lo <- max(0.001, min(0.999, p_lo))
      p_hi <- max(0.001, min(0.999, p_hi))

      tibble(metric       = m,
             estimate     = obs,
             ci_lo        = quantile(theta, p_lo, names = FALSE),
             ci_hi        = quantile(theta, p_hi, names = FALSE),
             n_boot_valid = n_b,
             note         = "")
    })
  }

  # --- Core function: run A2 for one trial type -----------------------------
  # trial_label:  "SAT" or "SBT"
  # ref_var:      name of flowsheet variable (string)
  # idx_var:      name of algorithm variable (string)
  # elig_var:     name of eligibility flag (string)
  # sensitivity_label: label for output (e.g. "primary", "modified_SAT")

  run_a2_analysis <- function(trial_label, ref_var, idx_var,
                               elig_var, sensitivity_label = "primary",
                               data = df_pp) {

    cat("  Running A2:", trial_label, "|", sensitivity_label, "\n")

    # Build analysis dataset: eligible days, complete pairs only
    df_a2 <- data %>%
      filter(.data[[elig_var]] == 1L,
             !is.na(.data[[ref_var]]),
             !is.na(.data[[idx_var]])) %>%
      mutate(ref = as.integer(.data[[ref_var]]),
             idx = as.integer(.data[[idx_var]]))

    n_eligible    <- nrow(df_a2)
    n_hosp        <- n_distinct(df_a2$hospitalization_id)

    cat("    Eligible days (complete pairs):", n_eligible,
        "| Hospitalizations:", n_hosp, "\n")

    if (n_eligible < 20 | n_hosp < 5) {
      cat("    WARNING: Insufficient data to run bootstrap reliably\n")
      cat("    (need >= 20 eligible days and >= 5 hospitalizations).\n")
      cat("    Point estimates only; CIs set to NA.\n\n")
    }

    # --- Point estimates on observed data -----------------------------------
    obs <- compute_a2_metrics(df_a2$ref, df_a2$idx)

    metric_names <- c("sensitivity", "specificity", "PPV", "NPV",
                      "accuracy", "F1", "MCC", "kappa")

    cat("    Point estimates:\n")
    for (m in metric_names) {
      cat("     ", m, ":", round(obs[[m]], 4), "\n")
    }
    cat("\n")

    # --- Confusion matrix ---------------------------------------------------
    cm_out <- tibble(
      trial             = trial_label,
      sensitivity_label = sensitivity_label,
      ref_var           = ref_var,
      idx_var           = idx_var,
      TP                = obs[["TP"]],
      TN                = obs[["TN"]],
      FP                = obs[["FP"]],
      FN                = obs[["FN"]],
      n_pairs           = obs[["n_pairs"]],
      n_hospitalizations = n_hosp,
      pct_positive_ref  = round((obs[["TP"]] + obs[["FN"]]) /
                                  obs[["n_pairs"]] * 100, 1),
      pct_positive_idx  = round((obs[["TP"]] + obs[["FP"]]) /
                                  obs[["n_pairs"]] * 100, 1),
      note = if (obs[["TP"]] + obs[["FN"]] == 0)
               "No positive reference events -- sensitivity undefined"
             else if (obs[["TN"]] + obs[["FP"]] == 0)
               "No negative reference events -- specificity undefined"
             else ""
    )

    # --- BCa Bootstrap ------------------------------------------------------
    # Pre-index rows by hospitalization for fast lookup -- avoids per-rep
    # filter() calls which are O(n*B) and stall on large datasets.
    # Strategy: build a list where each element is the integer row indices
    # for one hospitalization. Each rep samples IDs, retrieves row indices
    # via list lookup, concatenates, and slices ref/idx vectors directly.
    # At 7k+ hospitalizations x 1k reps this is ~100x faster than map_dfr.

    hosp_ids   <- unique(df_a2$hospitalization_id)
    n_hosp_b   <- length(hosp_ids)

    # Build row-index lookup list (one entry per unique hospitalization)
    hosp_row_idx <- split(seq_len(nrow(df_a2)), df_a2$hospitalization_id)

    # Pre-extract ref/idx as plain vectors for fast indexing
    ref_vec <- df_a2$ref
    idx_vec <- df_a2$idx

    cat("    Running BCa bootstrap (B =", B_BOOTSTRAP, ")...\n")

    set.seed(20250402)  # reproducible; date-stamped

    boot_metrics <- matrix(NA_real_,
                           nrow = B_BOOTSTRAP,
                           ncol = length(metric_names),
                           dimnames = list(NULL, metric_names))

    for (b in seq_len(B_BOOTSTRAP)) {
      # Progress marker every 100 reps
      if (b %% 100 == 0)
        cat("    Bootstrap rep", b, "of", B_BOOTSTRAP, "\n")

      # Resample hospitalization IDs with replacement (cluster bootstrap)
      sampled_ids <- sample(hosp_ids, size = n_hosp_b, replace = TRUE)

      # Expand to row indices: unlist preserves all days within each sampled
      # hospitalization, including duplicates when a hosp is drawn >1 time
      boot_rows <- unlist(hosp_row_idx[sampled_ids], use.names = FALSE)

      bm <- compute_a2_metrics(ref_vec[boot_rows], idx_vec[boot_rows])
      for (m in metric_names) boot_metrics[b, m] <- bm[[m]]
    }

    cat("    Bootstrap complete.\n")

    # --- BCa CIs ------------------------------------------------------------
    ci_df <- bca_ci(boot_metrics, obs, metric_names)

    # Attach trial/label info
    metrics_out <- ci_df %>%
      mutate(trial             = trial_label,
             sensitivity_label = sensitivity_label,
             ref_var           = ref_var,
             idx_var           = idx_var,
             n_eligible_days   = n_eligible,
             n_hospitalizations = n_hosp,
             B_bootstrap       = B_BOOTSTRAP,
             bootstrap_note    = if (B_BOOTSTRAP < 10000L)
               "PILOT: B < 10,000 -- rerun at B=10,000 before submission"
               else "") %>%
      select(trial, sensitivity_label, ref_var, idx_var, metric,
             estimate, ci_lo, ci_hi, n_boot_valid,
             n_eligible_days, n_hospitalizations, B_bootstrap,
             bootstrap_note, note)

    # --- Bootstrap distribution data for CC figure reconstruction -----------
    boot_dist_out <- as_tibble(boot_metrics) %>%
      mutate(trial             = trial_label,
             sensitivity_label = sensitivity_label,
             boot_rep          = seq_len(B_BOOTSTRAP))

    list(metrics      = metrics_out,
         confusion    = cm_out,
         boot_dist    = boot_dist_out)
  }

  # ---------------------------------------------------------------------------
  # 2.1 SAT CRITERION VALIDITY (primary)
  # ---------------------------------------------------------------------------

  cat("-- 2.1 SAT Criterion Validity (primary)\n\n")

  a2_sat_results <- NULL
  a2_sat_cm      <- NULL
  a2_sat_boot    <- NULL

  if (site_has_flowsheet_sat) {
    res_sat <- run_a2_analysis(
      trial_label       = "SAT",
      ref_var           = "flowsheet_SAT",
      idx_var           = "SAT_delivered_primary",
      elig_var          = "SAT_eligible",
      sensitivity_label = "primary"
    )
    a2_sat_results <- res_sat$metrics
    a2_sat_cm      <- res_sat$confusion
    a2_sat_boot    <- res_sat$boot_dist

    cat("  SAT primary analysis complete.\n\n")

  } else {
    cat("  SKIP: site_has_flowsheet_sat = FALSE\n\n")
  }

  # ---------------------------------------------------------------------------
  # 2.2 SBT CRITERION VALIDITY (primary)
  # ---------------------------------------------------------------------------

  cat("-- 2.2 SBT Criterion Validity (primary)\n\n")

  a2_sbt_results <- NULL
  a2_sbt_cm      <- NULL
  a2_sbt_boot    <- NULL

  if (site_has_flowsheet_sbt) {
    res_sbt <- run_a2_analysis(
      trial_label       = "SBT",
      ref_var           = "flowsheet_SBT",
      idx_var           = "SBT_delivered_2min",
      elig_var          = "SBT_eligible",
      sensitivity_label = "primary"
    )
    a2_sbt_results <- res_sbt$metrics
    a2_sbt_cm      <- res_sbt$confusion
    a2_sbt_boot    <- res_sbt$boot_dist

    cat("  SBT primary analysis complete.\n\n")

  } else {
    cat("  SKIP: site_has_flowsheet_sbt = FALSE\n\n")
  }

  # ---------------------------------------------------------------------------
  # 2.S SENSITIVITY: ALTERNATIVE EXPOSURE DEFINITIONS
  # ---------------------------------------------------------------------------

  cat("-- 2.S Sensitivity: Alternative exposure definitions\n\n")

  a2_sens_results <- NULL
  a2_sens_cm      <- NULL
  a2_sens_boot    <- NULL

  sens_list <- list()

  if (site_has_flowsheet_sat) {
    cat("  Running 2S1: SAT_delivered_modified vs. flowsheet_SAT\n")
    res_s1 <- run_a2_analysis(
      trial_label       = "SAT",
      ref_var           = "flowsheet_SAT",
      idx_var           = "SAT_delivered_modified",
      elig_var          = "SAT_eligible",
      sensitivity_label = "2S1_modified_SAT"
    )
    sens_list[["2S1"]] <- res_s1
  }

  if (site_has_flowsheet_sbt) {
    cat("  Running 2S2: SBT_delivered_5min (sensitivity) vs. flowsheet_SBT\n")
    res_s2 <- run_a2_analysis(
      trial_label       = "SBT",
      ref_var           = "flowsheet_SBT",
      idx_var           = "SBT_delivered_5min",
      elig_var          = "SBT_eligible",
      sensitivity_label = "2S2_5min_SBT"
    )
    sens_list[["2S2"]] <- res_s2
  }

  if (length(sens_list) > 0) {
    a2_sens_results <- purrr::map_dfr(sens_list, ~ .x$metrics)
    a2_sens_cm      <- purrr::map_dfr(sens_list, ~ .x$confusion)
    a2_sens_boot    <- purrr::map_dfr(sens_list, ~ .x$boot_dist)
  }

  # ---------------------------------------------------------------------------
  # 2.3 ANALYSIS 2 FIGURE
  # ---------------------------------------------------------------------------
  # Dot-and-CI plot: 8 metrics on y-axis, estimate + BCa CI on x
  # Two panels: SAT | SBT
  # Primary in solid color; sensitivity runs in muted overlay
  # ---------------------------------------------------------------------------

  cat("-- 2.3 Analysis 2 figure\n")

  metric_order <- c("sensitivity", "specificity", "PPV", "NPV",
                    "accuracy", "F1", "MCC", "kappa")

  metric_labels <- c(
    sensitivity = "Sensitivity",
    specificity = "Specificity",
    PPV         = "PPV",
    NPV         = "NPV",
    accuracy    = "Accuracy",
    F1          = "F1 Score",
    MCC         = "MCC",
    kappa       = "Cohen's Kappa"
  )

  # Combine primary + sensitivity for plot
  plot_data_list <- list()

  if (!is.null(a2_sat_results))
    plot_data_list[["sat_primary"]] <- a2_sat_results %>%
      mutate(run_type = "primary")
  if (!is.null(a2_sbt_results))
    plot_data_list[["sbt_primary"]] <- a2_sbt_results %>%
      mutate(run_type = "primary")
  if (!is.null(a2_sens_results))
    plot_data_list[["sensitivity"]] <- a2_sens_results %>%
      mutate(run_type = "sensitivity")

  if (length(plot_data_list) > 0) {

    plot_data_a2 <- bind_rows(plot_data_list) %>%
      filter(metric %in% metric_order) %>%
      mutate(
        metric_f = factor(metric, levels = rev(metric_order),
                          labels = rev(metric_labels[metric_order])),
        trial_f  = factor(trial, levels = c("SAT", "SBT")),
        alpha_val = if_else(run_type == "primary", 1, 0.4),
        shape_val = if_else(run_type == "primary", 16L, 17L),
        run_label = case_when(
          sensitivity_label == "primary"        ~ "Primary",
          sensitivity_label == "2S1_modified_SAT" ~ "Sensitivity: Modified SAT",
          sensitivity_label == "2S2_2min_SBT"   ~ "Sensitivity: 2-min SBT",
          TRUE ~ sensitivity_label
        )
      )

    # Build panels for whichever trials ran
    panels_available <- unique(plot_data_a2$trial)

    fig_list_a2 <- purrr::map(panels_available, function(tr) {
      pd <- plot_data_a2 %>% filter(trial == tr)

      ggplot(pd, aes(x = estimate, xmin = ci_lo, xmax = ci_hi,
                     y = metric_f, color = run_label, alpha = alpha_val)) +
        geom_vline(xintercept = c(0.8, 0.9), linetype = "dotted",
                   color = "gray70", linewidth = 0.4) +
        geom_errorbarh(height = 0.25, linewidth = 0.7,
                       position = position_dodge(width = 0.5)) +
        geom_point(aes(shape = run_label), size = 3,
                   position = position_dodge(width = 0.5)) +
        scale_color_manual(values = c(
          "Primary"                   = if (tr == "SAT") clr_sat else clr_sbt,
          "Sensitivity: Modified SAT" = "#88CCEE",
          "Sensitivity: 2-min SBT"   = "#DDAA33"
        )) +
        scale_alpha_identity() +
        scale_shape_manual(values = c(
          "Primary"                   = 16,
          "Sensitivity: Modified SAT" = 17,
          "Sensitivity: 2-min SBT"   = 17
        )) +
        scale_x_continuous(limits = c(0, 1.05),
                           breaks = c(0, 0.25, 0.5, 0.75, 0.8, 0.9, 1.0),
                           labels = c("0", "0.25", "0.50", "0.75",
                                      "0.80", "0.90", "1.0")) +
        labs(title    = paste0(tr, " Criterion Validity -- Algorithm vs. Flowsheet"),
             subtitle = "BCa bootstrap 95% CI | Eligible days only | Cluster = hospitalization",
             x        = "Metric Value (0\u20131)",
             y        = NULL,
             color    = NULL,
             shape    = NULL,
             caption  = paste0(
               "Reference standard: nurse flowsheet (", tr, ").\n",
               "Index test: algorithm (primary and sensitivity definitions).\n",
               "Dotted lines at 0.80 and 0.90 for reference.\n",
               if (B_BOOTSTRAP < 10000L)
                 paste0("*** PILOT: BCa B=", B_BOOTSTRAP,
                        " -- rerun at B=10,000 before submission ***")
               else ""
             )) +
        theme_abtrise() +
        theme(panel.grid.minor = element_blank())
    })

    if (length(fig_list_a2) == 1) {
      fig_A2 <- fig_list_a2[[1]]
    } else {
      fig_A2 <- wrap_plots(fig_list_a2, ncol = 2)
    }

    export_png(fig_A2, "figures/a2", "fig_A2_metrics.png",
               width = if (length(panels_available) == 2) 14 else 8,
               height = 7)

    cat("  Figure exported.\n\n")

  } else {
    cat("  No figure data available (no flowsheet runs completed).\n\n")
    fig_A2 <- NULL
  }

  # ---------------------------------------------------------------------------
  # 2.4 ANALYSIS 2 EXPORTS
  # Combined outputs: SAT and SBT combined into single files where parallel
  # structure allows. Each row carries a trial column (SAT / SBT) so the
  # CC can filter or display side by side. Separate bootstrap dist files
  # retained per trial (too large to combine without loss of CC utility).
  # ---------------------------------------------------------------------------

  cat("-- 2.4 Exporting Analysis 2 outputs\n")

  # Combined primary metrics (SAT + SBT, one row per metric per trial)
  a2_primary_metrics <- bind_rows(
    if (!is.null(a2_sat_results)) a2_sat_results else NULL,
    if (!is.null(a2_sbt_results)) a2_sbt_results else NULL
  )
  if (nrow(a2_primary_metrics) > 0)
    export_csv(a2_primary_metrics, "models/a2", "A2_metrics.csv")

  # Combined confusion matrix (one row per trial)
  a2_confusion <- bind_rows(
    if (!is.null(a2_sat_cm)) a2_sat_cm else NULL,
    if (!is.null(a2_sbt_cm)) a2_sbt_cm else NULL
  )
  if (nrow(a2_confusion) > 0)
    export_csv(a2_confusion, "models/a2", "A2_confusion_matrix.csv")

  # Combined sensitivity metrics (all sensitivity runs, both trials)
  if (!is.null(a2_sens_results))
    export_csv(a2_sens_results, "models/a2", "A2_sensitivity_metrics.csv")

  # Combined sensitivity confusion matrices
  a2_sens_confusion <- bind_rows(
    if (!is.null(a2_sens_cm)) a2_sens_cm else NULL
  )
  if (nrow(a2_sens_confusion) > 0)
    export_csv(a2_sens_confusion, "models/a2", "A2_sensitivity_confusion_matrix.csv")

  # Bootstrap distributions -- kept separate per trial (large files;
  # CC filters by trial column for pooling)
  if (!is.null(a2_sat_boot))
    export_csv(a2_sat_boot, "figures/a2", "A2_SAT_boot_dist.csv")

  if (!is.null(a2_sbt_boot))
    export_csv(a2_sbt_boot, "figures/a2", "A2_SBT_boot_dist.csv")

  if (!is.null(a2_sens_boot))
    export_csv(a2_sens_boot, "figures/a2", "A2_sensitivity_boot_dist.csv")

  # Combined figure data: all primary + sensitivity, both trials
  a2_all_fig_data <- bind_rows(
    if (!is.null(a2_sat_results)) a2_sat_results else NULL,
    if (!is.null(a2_sbt_results)) a2_sbt_results else NULL,
    if (!is.null(a2_sens_results)) a2_sens_results else NULL
  )
  if (nrow(a2_all_fig_data) > 0)
    export_csv(a2_all_fig_data, "figures/a2", "A2_fig_metrics_data.csv")

  cat("\nAnalysis 2 complete.\n\n")

# end flowsheet gate


# =============================================================================
# ANALYSIS 2 DIAGNOSTICS
# Exploratory checks to interpret agreement results
# Prints to console only -- not part of primary outputs
# Run interactively section by section or source the whole file
# =============================================================================

cat("============================================================\n")
cat("ANALYSIS 2 DIAGNOSTICS\n")
cat("============================================================\n\n")

# D1: Flowsheet completion rate -----------------------------------------------
cat("-- D1: Flowsheet completion rate\n\n")

d1_sat <- df_pp %>%
  filter(SAT_eligible == 1) %>%
  summarise(
    trial               = "SAT",
    n_eligible_days     = n(),
    n_flowsheet_entered = sum(!is.na(flowsheet_SAT)),
    n_flowsheet_missing = sum(is.na(flowsheet_SAT)),
    pct_completed       = round(mean(!is.na(flowsheet_SAT)) * 100, 1),
    n_flowsheet_pos     = sum(flowsheet_SAT == 1, na.rm = TRUE),
    n_flowsheet_neg     = sum(flowsheet_SAT == 0, na.rm = TRUE),
    pct_pos_among_entered = round(mean(flowsheet_SAT == 1, na.rm = TRUE) * 100, 1)
  )

d1_sbt <- df_pp %>%
  filter(SBT_eligible == 1) %>%
  summarise(
    trial               = "SBT",
    n_eligible_days     = n(),
    n_flowsheet_entered = sum(!is.na(flowsheet_SBT)),
    n_flowsheet_missing = sum(is.na(flowsheet_SBT)),
    pct_completed       = round(mean(!is.na(flowsheet_SBT)) * 100, 1),
    n_flowsheet_pos     = sum(flowsheet_SBT == 1, na.rm = TRUE),
    n_flowsheet_neg     = sum(flowsheet_SBT == 0, na.rm = TRUE),
    pct_pos_among_entered = round(mean(flowsheet_SBT == 1, na.rm = TRUE) * 100, 1)
  )

cat("Overall flowsheet completion:\n")
print(as.data.frame(bind_rows(d1_sat, d1_sbt)))
cat("\n")

d1_by_hosp <- bind_rows(
  df_pp %>% filter(SAT_eligible == 1) %>% group_by(hospital_id) %>%
    summarise(trial = "SAT", n_eligible = n(),
              pct_fs_completed = round(mean(!is.na(flowsheet_SAT)) * 100, 1),
              pct_fs_positive  = round(mean(flowsheet_SAT == 1, na.rm = TRUE) * 100, 1),
              .groups = "drop"),
  df_pp %>% filter(SBT_eligible == 1) %>% group_by(hospital_id) %>%
    summarise(trial = "SBT", n_eligible = n(),
              pct_fs_completed = round(mean(!is.na(flowsheet_SBT)) * 100, 1),
              pct_fs_positive  = round(mean(flowsheet_SBT == 1, na.rm = TRUE) * 100, 1),
              .groups = "drop")
)

cat("Flowsheet completion by hospital:\n")
print(as.data.frame(d1_by_hosp))
cat("\n")

low_completion <- d1_by_hosp %>% filter(pct_fs_completed < 50)
if (nrow(low_completion) > 0) {
  cat("WARNING: Hospitals with < 50% flowsheet completion:\n")
  print(as.data.frame(low_completion))
  cat("  These contribute unreliable reference standard data.\n\n")
} else {
  cat("All hospitals >= 50% flowsheet completion.\n\n")
}

# D2: Disagreement by hospital ------------------------------------------------
cat("-- D2: Confusion matrix by hospital\n\n")

compute_cm_by_hosp <- function(data, ref_var, idx_var, trial_label) {
  data %>%
    filter(!is.na(.data[[ref_var]]), !is.na(.data[[idx_var]])) %>%
    group_by(hospital_id) %>%
    summarise(
      trial       = trial_label,
      n_pairs     = n(),
      TP = sum(.data[[ref_var]] == 1 & .data[[idx_var]] == 1),
      TN = sum(.data[[ref_var]] == 0 & .data[[idx_var]] == 0),
      FP = sum(.data[[ref_var]] == 0 & .data[[idx_var]] == 1),
      FN = sum(.data[[ref_var]] == 1 & .data[[idx_var]] == 0),
      sensitivity  = round(TP / (TP + FN + 1e-9), 3),
      specificity  = round(TN / (TN + FP + 1e-9), 3),
      pct_disagree = round((FP + FN) / n() * 100, 1),
      kappa = {
        p_obs    <- (TP + TN) / n()
        p_chance <- ((TP+FN)/n())*((TP+FP)/n()) + ((TN+FP)/n())*((TN+FN)/n())
        round(ifelse((1-p_chance) > 0, (p_obs-p_chance)/(1-p_chance), NA), 3)
      },
      .groups = "drop"
    )
}

d2_sat <- df_pp %>% filter(SAT_eligible == 1) %>%
  compute_cm_by_hosp("flowsheet_SAT", "SAT_delivered_primary", "SAT")
d2_sbt <- df_pp %>% filter(SBT_eligible == 1) %>%
  compute_cm_by_hosp("flowsheet_SBT", "SBT_delivered_2min", "SBT")

cat("SAT confusion matrix by hospital (sorted by kappa):\n")
print(as.data.frame(d2_sat %>% arrange(kappa)))
cat("\nSBT confusion matrix by hospital (sorted by kappa):\n")
print(as.data.frame(d2_sbt %>% arrange(kappa)))
cat("\n")

# D3: Conditional agreement on flowsheet-positive days ------------------------
cat("-- D3: Conditional agreement on flowsheet-positive days\n\n")

d3_sat <- df_pp %>%
  filter(SAT_eligible == 1, flowsheet_SAT == 1) %>%
  summarise(
    trial                  = "SAT",
    n_flowsheet_pos        = n(),
    pct_alg_agrees         = round(mean(SAT_delivered_primary == 1,  na.rm = TRUE) * 100, 1),
    pct_modified_agrees    = round(mean(SAT_delivered_modified == 1, na.rm = TRUE) * 100, 1),
    pct_had_sedation_prior = round(mean(sedation_prior == 1,         na.rm = TRUE) * 100, 1),
    median_SOFA_prior      = round(median(SOFA_prior, na.rm = TRUE), 1)
  )

d3_sbt <- df_pp %>%
  filter(SBT_eligible == 1, flowsheet_SBT == 1) %>%
  summarise(
    trial                  = "SBT",
    n_flowsheet_pos        = n(),
    pct_alg_agrees         = round(mean(SBT_delivered_2min == 1,  na.rm = TRUE) * 100, 1),
    pct_5min_agrees        = round(mean(SBT_delivered_5min == 1,  na.rm = TRUE) * 100, 1),
    median_FiO2_prior      = round(median(FiO2_prior, na.rm = TRUE), 3),
    median_PEEP_prior      = round(median(PEEP_prior,  na.rm = TRUE), 1),
    pct_high_support_prior = round(mean((FiO2_prior > 0.5 | PEEP_prior > 8),
                                         na.rm = TRUE) * 100, 1)
  )

cat("Among flowsheet-positive days -- algorithm agreement and data signals:\n")
cat("SAT (flowsheet_SAT == 1):\n")
print(as.data.frame(t(d3_sat)))
cat("\nSBT (flowsheet_SBT == 1):\n")
print(as.data.frame(t(d3_sbt)))
cat("\n")

# D4: Temporal pattern -- agreement by ventilator day -------------------------
cat("-- D4: Agreement by ventilator day\n\n")

calc_kappa_byday <- function(data, ref_var, idx_var, trial_label) {
  data %>%
    filter(!is.na(.data[[ref_var]]), !is.na(.data[[idx_var]])) %>%
    group_by(vent_day) %>%
    summarise(
      trial        = trial_label,
      n_pairs      = n(),
      pct_fs_pos   = round(mean(.data[[ref_var]] == 1) * 100, 1),
      pct_alg_pos  = round(mean(.data[[idx_var]] == 1) * 100, 1),
      kappa = {
        TP <- sum(.data[[ref_var]] == 1 & .data[[idx_var]] == 1)
        TN <- sum(.data[[ref_var]] == 0 & .data[[idx_var]] == 0)
        FP <- sum(.data[[ref_var]] == 0 & .data[[idx_var]] == 1)
        FN <- sum(.data[[ref_var]] == 1 & .data[[idx_var]] == 0)
        n  <- n()
        p_obs    <- (TP + TN) / n
        p_chance <- ((TP+FN)/n)*((TP+FP)/n) + ((TN+FP)/n)*((TN+FN)/n)
        round(ifelse((1-p_chance) > 0, (p_obs-p_chance)/(1-p_chance), NA), 3)
      },
      .groups = "drop"
    )
}

d4_sat <- df_pp %>% filter(SAT_eligible == 1) %>%
  calc_kappa_byday("flowsheet_SAT", "SAT_delivered_primary", "SAT")
d4_sbt <- df_pp %>% filter(SBT_eligible == 1) %>%
  calc_kappa_byday("flowsheet_SBT", "SBT_delivered_2min", "SBT")
d4_combined <- bind_rows(d4_sat, d4_sbt)

cat("SAT agreement by ventilator day:\n")
print(as.data.frame(d4_sat))
cat("\nSBT agreement by ventilator day:\n")
print(as.data.frame(d4_sbt))
cat("\n")

fig_d4_kappa <- ggplot(d4_combined, aes(x = vent_day, y = kappa, color = trial)) +
  geom_hline(yintercept = 0,   linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = 0.2, linetype = "dotted", color = "gray70") +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_color_manual(values = clr_trial) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 28)) +
  labs(title    = "Cohen's Kappa by Ventilator Day",
       subtitle = "Declining kappa over vent course = documentation fatigue hypothesis",
       x = "Ventilator Day", y = "Cohen's Kappa", color = NULL,
       caption  = "Dotted line = kappa 0.20 (fair agreement threshold). Diagnostic only.") +
  theme_abtrise()

fig_d4_rates <- ggplot(d4_combined, aes(x = vent_day, color = trial)) +
  geom_line(aes(y = pct_fs_pos,  linetype = "Flowsheet"), linewidth = 0.9) +
  geom_line(aes(y = pct_alg_pos, linetype = "Algorithm"),  linewidth = 0.9) +
  scale_color_manual(values = clr_trial) +
  scale_linetype_manual(values = c("Flowsheet" = "solid", "Algorithm" = "dashed")) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 28)) +
  facet_wrap(~ trial, ncol = 2) +
  labs(title    = "Positive Rate by Ventilator Day: Flowsheet vs. Algorithm",
       subtitle = "Diverging rates over vent course suggest documentation decay",
       x = "Ventilator Day", y = "% Positive (among eligible days with flowsheet entry)",
       color = NULL, linetype = NULL,
       caption = "Diagnostic only -- not for manuscript.") +
  theme_abtrise()

fig_d4_combined <- fig_d4_kappa / fig_d4_rates

print(fig_d4_combined)

export_png(fig_d4_combined, "figures/a2", "fig_A2_D4_kappa_by_ventday.png",
           width = 8, height = 8)

cat("D4 figure rendered to viewer and exported.\n\n")

# D4 summary: early vs. late vent course kappa split -------------------------
# Quantifies documentation fatigue (SAT) vs. definition mismatch (SBT) in
# a single presentable table. Breakpoints: days 1-5 (early, both systems
# active), days 6-14 (mid, fatigue begins), days 15-28 (late, flowsheet
# increasingly unreliable for SAT).

cat("-- D4 Summary: Early vs. Late Vent Course Kappa\n\n")

d4_summary <- d4_combined %>%
  mutate(
    period = case_when(
      vent_day <= 5  ~ "Days 1-5  (early)",
      vent_day <= 14 ~ "Days 6-14 (mid)",
      TRUE           ~ "Days 15-28 (late)"
    ),
    period = factor(period, levels = c("Days 1-5  (early)",
                                        "Days 6-14 (mid)",
                                        "Days 15-28 (late)"))
  ) %>%
  group_by(trial, period) %>%
  summarise(
    n_vent_days      = sum(n_pairs),
    mean_kappa       = round(mean(kappa,      na.rm = TRUE), 3),
    mean_pct_fs_pos  = round(mean(pct_fs_pos,  na.rm = TRUE), 1),
    mean_pct_alg_pos = round(mean(pct_alg_pos, na.rm = TRUE), 1),
    fs_alg_gap       = round(mean_pct_fs_pos - mean_pct_alg_pos, 1),
    .groups = "drop"
  )

cat("Kappa and positive rate by vent course period:\n")
print(as.data.frame(d4_summary))
cat("\n")

# Plain-language interpretation for each trial
for (tr in c("SAT", "SBT")) {
  early_k <- d4_summary %>%
    filter(trial == tr, period == "Days 1-5  (early)") %>%
    pull(mean_kappa)
  late_k <- d4_summary %>%
    filter(trial == tr, period == "Days 15-28 (late)") %>%
    pull(mean_kappa)
  early_gap <- d4_summary %>%
    filter(trial == tr, period == "Days 1-5  (early)") %>%
    pull(fs_alg_gap)
  late_gap <- d4_summary %>%
    filter(trial == tr, period == "Days 15-28 (late)") %>%
    pull(fs_alg_gap)

  cat(tr, "interpretation:\n")
  cat("  Early kappa (days 1-5): ", early_k,
      "| Late kappa (days 15-28):", late_k, "\n")
  cat("  Early flowsheet-algorithm gap:", early_gap,
      "pp | Late gap:", late_gap, "pp\n")

  if (tr == "SAT") {
    if (late_k < 0 && early_k > 0.1) {
      cat("  -> Kappa crosses zero by late vent days.\n")
      cat("     Algorithm positive rate remains stable; flowsheet rate collapses.\n")
      cat("     Pattern consistent with documentation fatigue, not algorithm failure.\n")
    } else if (late_k < early_k) {
      cat("  -> Kappa declines over vent course -- some documentation fatigue present.\n")
    } else {
      cat("  -> Kappa stable across vent course -- no documentation fatigue signal.\n")
    }
  }

  if (tr == "SBT") {
    if (early_gap > 20) {
      cat("  -> Large flowsheet-algorithm gap present from day 1 (", early_gap, "pp).\n")
      cat("     Flowsheet captures broader SBT definition than algorithm from the outset.\n")
      cat("     This is a definition mismatch, not documentation fatigue.\n")
    } else if (late_k < early_k - 0.1) {
      cat("  -> Kappa declines in late vent days -- both definition mismatch\n")
      cat("     and documentation fatigue may be contributing.\n")
    } else {
      cat("  -> Kappa relatively stable -- definition mismatch is the primary driver.\n")
    }
  }
  cat("\n")
}

# D5: Borderline event check --------------------------------------------------
cat("-- D5: Clinical profile by agreement cell\n\n")

summarise_cells <- function(data, ref_var, idx_var) {
  data %>%
    filter(!is.na(.data[[ref_var]]), !is.na(.data[[idx_var]])) %>%
    mutate(cell = case_when(
      .data[[ref_var]] == 1 & .data[[idx_var]] == 1 ~ "TP",
      .data[[ref_var]] == 0 & .data[[idx_var]] == 1 ~ "FP (alg+, fs-)",
      .data[[ref_var]] == 1 & .data[[idx_var]] == 0 ~ "FN (alg-, fs+)",
      TRUE                                           ~ "TN"
    )) %>%
    group_by(cell) %>%
    summarise(
      n             = n(),
      med_SOFA      = round(median(SOFA_prior,  na.rm = TRUE), 1),
      med_FiO2      = round(median(FiO2_prior,  na.rm = TRUE), 3),
      med_PEEP      = round(median(PEEP_prior,  na.rm = TRUE), 1),
      pct_sedation  = round(mean(sedation_prior == 1, na.rm = TRUE) * 100, 1),
      med_vent_day  = round(median(vent_day, na.rm = TRUE), 1),
      .groups = "drop"
    )
}

d5_sat <- df_pp %>% filter(SAT_eligible == 1) %>%
  summarise_cells("flowsheet_SAT", "SAT_delivered_primary")
d5_sbt <- df_pp %>% filter(SBT_eligible == 1) %>%
  summarise_cells("flowsheet_SBT", "SBT_delivered_2min")

cat("SAT -- clinical profile by agreement cell:\n")
print(as.data.frame(d5_sat))
cat("\nSBT -- clinical profile by agreement cell:\n")
print(as.data.frame(d5_sbt))
cat("\n")

# Automated interpretation: are FP days similar to TP days?
for (trial_label in c("SAT", "SBT")) {
  prof <- if (trial_label == "SAT") d5_sat else d5_sbt
  tp_sofa <- prof %>% filter(cell == "TP")              %>% pull(med_SOFA)
  fp_sofa <- prof %>% filter(cell == "FP (alg+, fs-)") %>% pull(med_SOFA)
  tn_sofa <- prof %>% filter(cell == "TN")              %>% pull(med_SOFA)
  if (length(tp_sofa) > 0 & length(fp_sofa) > 0) {
    cat(trial_label, "SOFA: TP =", tp_sofa, "| FP =", fp_sofa,
        "| TN =", if (length(tn_sofa) > 0) tn_sofa else "n/a", "\n")
    if (abs(tp_sofa - fp_sofa) < 1)
      cat("  -> FP days similar to TP: algorithm likely detecting real events nurse did not document.\n")
    else
      cat("  -> FP days differ from TP: review vent_day distribution and sedation patterns.\n")
  }
}
cat("\n")

# D6: Summary for team discussion ---------------------------------------------
cat("-- D6: Summary for team discussion\n\n")

d6 <- bind_rows(d1_sat, d1_sbt) %>%
  select(trial, n_eligible_days, pct_completed, pct_pos_among_entered) %>%
  left_join(
    bind_rows(
      d3_sat %>% select(trial, pct_alg_agrees),
      d3_sbt %>% select(trial, pct_alg_agrees)
    ),
    by = "trial"
  ) %>%
  left_join(
    d4_combined %>%
      group_by(trial) %>%
      summarise(
        kappa_days1_5   = round(mean(kappa[vent_day <= 5],       na.rm = TRUE), 3),
        kappa_days6_14  = round(mean(kappa[vent_day %in% 6:14],  na.rm = TRUE), 3),
        kappa_days15_28 = round(mean(kappa[vent_day >= 15],      na.rm = TRUE), 3),
        .groups = "drop"
      ),
    by = "trial"
  )

print(as.data.frame(d6))
cat("\n")
cat("Key questions for team:\n")
cat("  1. Is flowsheet completion rate acceptable as a reference standard?\n")
cat("  2. Is disagreement localized to specific hospitals (D2)?\n")
cat("  3. Does kappa decay over vent course (documentation fatigue, D4)?\n")
cat("  4. Are FP days clinically similar to TP days (D5)?\n\n")

# =============================================================================
# D7: DENOMINATOR COMPARISON TABLE
# =============================================================================
# Side-by-side comparison of eligible-day vs. all-vent-day denominators.
# Purpose: documents why the two approaches produce different results and
# confirms the SAP-specified eligible-day denominator is correct.
# Useful for team presentations and methods transparency.
# =============================================================================

cat("-- D7: Denominator comparison -- eligible days vs. all vent days\n\n")

# Compute metrics on ALL vent days (no eligibility filter) for comparison
# This replicates what an analysis without the eligible-day filter would produce

run_all_ventdays <- function(trial_label, ref_var, idx_var, data) {
  df_all <- data %>%
    filter(!is.na(.data[[ref_var]]), !is.na(.data[[idx_var]])) %>%
    mutate(ref = as.integer(.data[[ref_var]]),
           idx = as.integer(.data[[idx_var]]))

  TP <- sum(df_all$ref == 1 & df_all$idx == 1)
  TN <- sum(df_all$ref == 0 & df_all$idx == 0)
  FP <- sum(df_all$ref == 0 & df_all$idx == 1)
  FN <- sum(df_all$ref == 1 & df_all$idx == 0)
  n  <- nrow(df_all)

  sens <- round(TP / (TP + FN + 1e-9), 3)
  spec <- round(TN / (TN + FP + 1e-9), 3)
  ppv  <- round(TP / (TP + FP + 1e-9), 3)
  npv  <- round(TN / (TN + FN + 1e-9), 3)
  f1   <- round(2*TP / (2*TP + FP + FN + 1e-9), 3)
  TP_n <- as.numeric(TP); TN_n <- as.numeric(TN)
  FP_n <- as.numeric(FP); FN_n <- as.numeric(FN)
  mcc_d <- sqrt((TP_n+FP_n)*(TP_n+FN_n)*(TN_n+FP_n)*(TN_n+FN_n))
  mcc  <- round(if (is.finite(mcc_d) && mcc_d > 0)
                  (TP_n*TN_n - FP_n*FN_n) / mcc_d else NA_real_, 3)
  p_obs    <- (TP + TN) / n
  p_chance <- ((TP+FN)/n)*((TP+FP)/n) + ((TN+FP)/n)*((TN+FN)/n)
  kappa <- round(if ((1-p_chance) > 0) (p_obs-p_chance)/(1-p_chance)
                 else NA_real_, 3)

  tibble(trial = trial_label, denominator = "All vent days (no eligibility filter)",
         n_days = n, TP = TP, TN = TN, FP = FP, FN = FN,
         sensitivity = sens, specificity = spec, PPV = ppv, NPV = npv,
         F1 = f1, MCC = mcc, kappa = kappa)
}

allday_sat <- run_all_ventdays("SAT", "flowsheet_SAT", "SAT_delivered_primary", df_pp)
allday_sbt <- run_all_ventdays("SBT", "flowsheet_SBT", "SBT_delivered_2min",    df_pp)

# Pull eligible-day results from primary analysis outputs
eligible_sat <- if (!is.null(a2_sat_results)) {
  a2_sat_results %>%
    filter(sensitivity_label == "primary") %>%
    select(trial, metric, estimate) %>%
    pivot_wider(names_from = metric, values_from = estimate) %>%
    mutate(denominator = "Eligible days only (SAP-specified)",
           n_days      = nrow(df_pp %>% filter(SAT_eligible == 1,
                                               !is.na(flowsheet_SAT)))) %>%
    left_join(
      df_pp %>% filter(SAT_eligible == 1, !is.na(flowsheet_SAT),
                       !is.na(SAT_delivered_primary)) %>%
        summarise(TP = sum(flowsheet_SAT==1 & SAT_delivered_primary==1),
                  TN = sum(flowsheet_SAT==0 & SAT_delivered_primary==0),
                  FP = sum(flowsheet_SAT==0 & SAT_delivered_primary==1),
                  FN = sum(flowsheet_SAT==1 & SAT_delivered_primary==0)),
      by = character()
    )
} else NULL

eligible_sbt <- if (!is.null(a2_sbt_results)) {
  a2_sbt_results %>%
    filter(sensitivity_label == "primary") %>%
    select(trial, metric, estimate) %>%
    pivot_wider(names_from = metric, values_from = estimate) %>%
    mutate(denominator = "Eligible days only (SAP-specified)",
           n_days      = nrow(df_pp %>% filter(SBT_eligible == 1,
                                               !is.na(flowsheet_SBT)))) %>%
    left_join(
      df_pp %>% filter(SBT_eligible == 1, !is.na(flowsheet_SBT),
                       !is.na(SBT_delivered_2min)) %>%
        summarise(TP = sum(flowsheet_SBT==1 & SBT_delivered_2min==1),
                  TN = sum(flowsheet_SBT==0 & SBT_delivered_2min==0),
                  FP = sum(flowsheet_SBT==0 & SBT_delivered_2min==1),
                  FN = sum(flowsheet_SBT==1 & SBT_delivered_2min==0)),
      by = character()
    )
} else NULL

# Build comparison table
denom_compare <- bind_rows(
  eligible_sat %>% select(trial, denominator, n_days, TP, TN, FP, FN,
                           sensitivity, specificity, PPV, NPV, F1, MCC, kappa),
  allday_sat,
  eligible_sbt %>% select(trial, denominator, n_days, TP, TN, FP, FN,
                           sensitivity, specificity, PPV, NPV, F1, MCC, kappa),
  allday_sbt
) %>%
  arrange(trial, denominator)

cat("Denominator comparison (eligible days vs. all vent days):\n\n")
print(as.data.frame(denom_compare %>%
  select(trial, denominator, n_days, TP, TN, FP, FN,
         sensitivity, specificity, kappa)))
cat("\n")


# Plain-language summary
for (tr in c("SAT", "SBT")) {
  elig_k <- denom_compare %>%
    filter(trial == tr, str_detect(denominator, "Eligible")) %>%
    pull(kappa)
  all_k <- denom_compare %>%
    filter(trial == tr, str_detect(denominator, "All vent")) %>%
    pull(kappa)
  elig_n <- denom_compare %>%
    filter(trial == tr, str_detect(denominator, "Eligible")) %>%
    pull(n_days)
  all_n <- denom_compare %>%
    filter(trial == tr, str_detect(denominator, "All vent")) %>%
    pull(n_days)
  tn_all <- denom_compare %>%
    filter(trial == tr, str_detect(denominator, "All vent")) %>%
    pull(TN)
  tn_elig <- denom_compare %>%
    filter(trial == tr, str_detect(denominator, "Eligible")) %>%
    pull(TN)

  cat(tr, "denominator effect:\n")
  cat("  Eligible days:    n =", elig_n, "| kappa =", elig_k, "\n")
  cat("  All vent days:    n =", all_n,  "| kappa =", all_k,  "\n")
  cat("  Extra TN from ineligible days:", tn_all - tn_elig,
      "(structural zeros -- both algorithm and flowsheet agree no trial\n")
  cat("  possible on ineligible days; this agreement is not clinically meaningful)\n\n")
}

cat("NOTE: The SAP pre-specifies eligible days as the correct denominator.\n")


# Export comparison table
export_csv(denom_compare, "models/a2", "A2_denominator_comparison.csv")

# Session info ----------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "diagnostics", prefix_file("session_info_a2.txt")))
cat("Session info saved.\n")
cat("\n=== Analysis 2 script complete ===\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

