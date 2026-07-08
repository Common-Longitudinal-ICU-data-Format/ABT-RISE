# =============================================================================
# ABT-RISE: Run-All Wrapper Script
#
# PURPOSE:
#   Runs all four site-level analysis scripts in sequence with a single
#   Source click. Each script sources ABTRISE_01_setup.R automatically,
#   so data is loaded fresh for each analysis.
#
# HOW TO RUN:
#   1. Edit clif_config.json at the repo root (site_name, abtrise_input_dir,
#      abtrise_output_dir). No edits needed in this file.
#   2. Open this file in RStudio
#   3. Click Source (or Ctrl+Shift+Enter)
#   4. All outputs will appear in output_to_share/ subfolders
#
# SCRIPTS RUN IN ORDER:
#   1. ABTRISE_01_setup.R       -- data load, diagnostics, Table 1
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
# ESTIMATED RUNTIME (single-hospital site):
#   ~5-15 min depending on data size. Multi-hospital sites or B=10000
#   will take longer, primarily due to the bootstrap in ABTRISE_02.
#
# COORDINATING CENTER: Rush University Medical Center
# =============================================================================

# =============================================================================
# SECTION 0: SITE CONFIGURATION
# =============================================================================
# Site config (site_id, data_dir, out_dir) is read from clif_config.json at
# the repo root via ABTRISE_config.R. Sites edit clif_config.json only --
# no edits needed in this file.

# `here` is needed to locate config.R below, which then defines ensure_packages()
# for the rest of the pipeline. Bootstrap it first (runs under --vanilla, so pass
# an explicit CRAN mirror).
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(here))
source(here::here("code", "ABTRISE_config.R"))

# ANSI color codes for highlighting failures in the terminal. Always emitted
# (run_all.sh pipes stdout through tee; modern terminals render the codes,
# while log viewers like `less -R` handle them transparently).
.red   <- "\033[1;31m"
.reset <- "\033[0m"

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
  "ABTRISE_01_setup_c.R",
  "ABTRISE_02_criterion_c.R",
  "ABTRISE_345_outcomes_c.R",
  "ABTRISE_06_benchmarking_c.R"
)

results <- list()

# =============================================================================
# SECTION 2: PRE-FLIGHT CHECK -- ALL SCRIPTS MUST BE PRESENT
# =============================================================================
# Check that all four analysis scripts exist before running anything.
# If any are missing, stop immediately with a clear message listing what
# is missing and where to get it. Do not proceed with partial outputs.

cat("Pre-flight check: verifying all required scripts are present...\n")

missing_scripts <- scripts[!file.exists(here::here("code", scripts))]

if (length(missing_scripts) > 0) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                 MISSING REQUIRED SCRIPTS                     ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")
  cat("The following scripts are missing from your project folder:\n\n")
  for (s in missing_scripts) {
    cat("  ✗  ", s, "\n")
  }
  cat("\nAll five files must be in the code/ folder before running:\n")
  cat("  - code/ABTRISE_run_all.R          (this file)\n")
  cat("  - code/ABTRISE_01_setup_c.R\n")
  cat("  - code/ABTRISE_02_criterion_c.R\n")
  cat("  - code/ABTRISE_345_outcomes_c.R\n")
  cat("  - code/ABTRISE_06_benchmarking_c.R\n")
  cat("\nDownload all files from the project repository and save them\n")
  cat("in the code/ folder before running ABTRISE_run_all.R.\n\n")
  stop("Run aborted: missing required scripts. See message above.", call. = FALSE)
}

cat("  All required scripts found.\n\n")

# =============================================================================
# SECTION 3: RUN SCRIPTS IN ORDER
# =============================================================================

for (script in scripts) {

  script_path <- here::here("code", script)

  cat("────────────────────────────────────────────────────────────\n")
  cat("▶  Running:", script, "\n")
  cat("   Started:", format(Sys.time(), "%H:%M:%S"), "\n")
  cat("────────────────────────────────────────────────────────────\n\n")

  t_start <- proc.time()

  # withCallingHandlers preserves the "muffleWarning" restart so warnings
  # raised inside the sourced script don't interrupt execution; tryCatch
  # catches hard errors and lets the next script continue.
  tryCatch(
    withCallingHandlers(
      {
        source(script_path, local = FALSE)
        elapsed <- round((proc.time() - t_start)["elapsed"] / 60, 1)
        cat("\n✓  COMPLETE:", script, "--", elapsed, "min\n\n")
        results[[script]] <- list(status  = "success",
                                  message = "Completed without error",
                                  duration_min = elapsed)
      },
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) {
      elapsed <- round((proc.time() - t_start)["elapsed"] / 60, 1)
      cat("\n", .red, "✗  ERROR in ", script, " after ", elapsed, " min:",
          .reset, "\n", sep = "")
      cat("   ", .red, conditionMessage(e), .reset, "\n\n", sep = "")
      cat("   Continuing to next script...\n\n")
      results[[script]] <<- list(status  = "error",
                                 message = conditionMessage(e),
                                 duration_min = elapsed)
    }
  )
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

    "?"
  )
  dur <- if (!is.na(res$duration_min)) paste0(res$duration_min, " min") else "n/a"
  row_open  <- if (identical(res$status, "error")) .red   else ""
  row_close <- if (identical(res$status, "error")) .reset else ""
  cat(sprintf("║  %s%s%s  %-38s  %7s  ║\n",
              row_open, icon, row_close,
              substr(script, 1, 38),
              dur))
}

cat("╠══════════════════════════════════════════════════════════════╣\n")
cat("║  Total runtime: ", formatC(paste0(total_min, " min"),
                                   width = 44, flag = "-"), "║\n")
cat("║  Finished:      ", formatC(format(run_end, "%Y-%m-%d %H:%M:%S"),
                                   width = 44, flag = "-"), "║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Flag any failures (entire warning block printed in red)
errors <- Filter(function(r) r$status == "error", results)
if (length(errors) > 0) {
  cat(.red, "⚠  WARNING: ", length(errors), " script(s) failed:", .reset,
      "\n", sep = "")
  for (nm in names(errors)) {
    cat("   ", .red, "- ", nm, ": ", errors[[nm]]$message, .reset, "\n",
        sep = "")
  }
  cat("\n   ", .red, "Check ", out_dir,
      " to confirm which files were produced.", .reset, "\n", sep = "")
  cat("   ", .red, "Contact the coordinating center if errors persist.",
      .reset, "\n\n", sep = "")
} else {
  cat("All scripts completed successfully.\n")
  cat("Outputs are in:", out_dir, "\n\n")
}
