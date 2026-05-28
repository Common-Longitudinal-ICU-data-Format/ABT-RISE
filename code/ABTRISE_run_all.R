# =============================================================================
# ABT-RISE: Run-All Wrapper Script
#
# PURPOSE:
#   Runs all four site-level analysis scripts in sequence with a single
#   Source click. Each script sources ABTRISE_00_setup.R automatically,
#   so data is loaded fresh for each analysis.
#
# HOW TO RUN:
#   1. Edit clif_config.json at repo root (site_name, abtrise_input_dir,
#      abtrise_output_dir).
#   2. Open this file in RStudio.
#   3. Click Source (or Ctrl+Shift+Enter).
#   4. All outputs will appear in the configured output directory.
#
# SCRIPTS RUN IN ORDER:
#   1. ABTRISE_00_setup.R       -- data load, diagnostics, Table 1
#   2. ABTRISE_02_criterion.R   -- criterion validity (A2)
#   3. ABTRISE_345_outcomes.R   -- construct validity (A3, A4, A5, SA)
#   4. ABTRISE_06_benchmarking.R -- hospital benchmarking (A6)
#
# ERROR HANDLING:
#   - Each script is wrapped in tryCatch so a failure in one script does
#     not prevent the others from running
#   - Errors are logged to the console with the script name and message
#   - Check outputs/ after running to confirm all expected files were created
#
# ESTIMATED RUNTIME (single-hospital site, pilot B=1000):
#   ~5-15 min depending on data size. Multi-hospital sites or B=10000
#   will take longer, primarily due to the bootstrap in ABTRISE_02.
#
# COORDINATING CENTER: Rush University Medical Center
# =============================================================================

# =============================================================================
# SECTION 0: SITE CONFIGURATION
# =============================================================================
# Sites edit clif_config.json at repo root, not this file.
# Helper loads site_id, data_dir, out_dir from that config.

suppressPackageStartupMessages(library(here))
source(here::here("code", "ABTRISE_config.R"))

# =============================================================================
# SECTION 1: SETUP
# =============================================================================

run_start <- Sys.time()

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║              ABT-RISE Site Analysis -- Run All               ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat("║  Site:   ", formatC(site_id, width = 50, flag = "-"), "║\n")
cat("║  Start:  ", formatC(format(run_start, "%Y-%m-%d %H:%M:%S"),
                            width = 50, flag = "-"), "║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

scripts <- c(
  "ABTRISE_00_setup_c.R",
  "ABTRISE_02_criterion_c.R",
  "ABTRISE_345_outcomes_c.R",
  "ABTRISE_06_benchmarking_c.R"
)

results <- list()

# =============================================================================
# SECTION 2: RUN SCRIPTS IN ORDER
# =============================================================================

for (script in scripts) {

  script_path <- here::here("code", script)

  cat("────────────────────────────────────────────────────────────\n")
  cat("▶  Running:", script, "\n")
  cat("   Started:", format(Sys.time(), "%H:%M:%S"), "\n")
  cat("────────────────────────────────────────────────────────────\n\n")

  if (!file.exists(script_path)) {
    msg <- paste0("File not found: ", script_path)
    cat("✗  SKIPPED: ", msg, "\n\n")
    results[[script]] <- list(status = "skipped", message = msg,
                              duration_min = NA)
    next
  }

  t_start <- proc.time()

  # withCallingHandlers catches warnings WITHOUT unwinding the call stack,
  # so inner withCallingHandlers + invokeRestart("muffleWarning") calls in
  # the analysis scripts remain valid. tryCatch alone would unwind the stack
  # and invalidate those restarts, causing the "no restart muffleWarning found"
  # error. The outer tryCatch here catches only hard errors.
  tryCatch(
    withCallingHandlers(
      source(script_path, local = FALSE),
      warning = function(w) {
        # Warnings are already captured within each analysis script.
        # Muffle here so they don't print twice at the run_all level.
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      elapsed <- round((proc.time() - t_start)["elapsed"] / 60, 1)
      cat("\n\u2717  ERROR in", script, "after", elapsed, "min:\n")
      cat("   ", conditionMessage(e), "\n\n")
      cat("   Continuing to next script...\n\n")
      results[[script]] <<- list(status  = "error",
                                 message = conditionMessage(e),
                                 duration_min = elapsed)
    }
  )

  # Log success if not already logged as error
  if (is.null(results[[script]])) {
    elapsed <- round((proc.time() - t_start)["elapsed"] / 60, 1)
    cat("\n\u2713  COMPLETE:", script, "--", elapsed, "min\n\n")
    results[[script]] <- list(status       = "success",
                              message      = "Completed without error",
                              duration_min = elapsed)
  }
}

# =============================================================================
# SECTION 3: RUN SUMMARY
# =============================================================================

run_end  <- Sys.time()
total_min <- round(as.numeric(difftime(run_end, run_start, units = "mins")), 1)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║                     RUN SUMMARY                              ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")

for (script in scripts) {
  res <- results[[script]]
  icon <- switch(res$status,
    success = "✓",
    error   = "✗",
    skipped = "—",
    "?"
  )
  dur <- if (!is.na(res$duration_min)) paste0(res$duration_min, " min") else "n/a"
  cat(sprintf("║  %s  %-38s  %7s  ║\n",
              icon,
              substr(script, 1, 38),
              dur))
}

cat("╠══════════════════════════════════════════════════════════════╣\n")
cat("║  Total runtime: ", formatC(paste0(total_min, " min"),
                                   width = 44, flag = "-"), "║\n")
cat("║  Finished:      ", formatC(format(run_end, "%Y-%m-%d %H:%M:%S"),
                                   width = 44, flag = "-"), "║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Flag any failures or skips
errors  <- Filter(function(r) r$status == "error",   results)
skipped <- Filter(function(r) r$status == "skipped", results)

if (length(errors) > 0) {
  cat("⚠  WARNING:", length(errors), "script(s) failed:\n")
  for (nm in names(errors)) {
    cat("   -", nm, ":", errors[[nm]]$message, "\n")
  }
}
if (length(skipped) > 0) {
  cat("⚠  WARNING:", length(skipped), "script(s) skipped:\n")
  for (nm in names(skipped)) {
    cat("   -", nm, ":", skipped[[nm]]$message, "\n")
  }
}

if (length(errors) == 0 && length(skipped) == 0) {
  cat("All scripts completed successfully.\n")
  cat("Outputs are in: ", out_dir, "\n\n", sep = "")
} else {
  cat("\n   Check ", out_dir, " to confirm which files were produced.\n", sep = "")
  cat("   Contact the coordinating center if errors persist.\n\n")
  # Non-zero exit so the shell wrapper (run_all.sh / run_all.ps1) flags it.
  quit(status = 1, save = "no")
}
