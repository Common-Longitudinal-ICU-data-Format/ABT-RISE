#!/usr/bin/env bash
# ABT-RISE pipeline runner (macOS/Linux)
#
# Runs 5 Python steps then the R orchestrator. Tees everything to
# run_log_YYYYMMDD_HHMM.log at the repo root. Skips R gracefully and
# prints RStudio fallback instructions if Rscript is not on PATH.
#
# Usage:  chmod +x run_all.sh && ./run_all.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

LOG_FILE="$PROJECT_ROOT/run_log_$(date +%Y%m%d_%H%M).log"

# Capture all stdout+stderr of the rest of the script into the log.
exec > >(tee -a "$LOG_FILE") 2>&1

echo "ABT-RISE pipeline"
echo "Started: $(date)"
echo "Log:     $LOG_FILE"
echo

# ── config check ──────────────────────────────────────────────────────────
if [ ! -f "$PROJECT_ROOT/clif_config.json" ]; then
  echo "ERROR: clif_config.json not found at repo root."
  echo "       Copy clif_config_template.json to clif_config.json and edit it."
  exit 1
fi

# ── tool checks ───────────────────────────────────────────────────────────
if ! command -v uv >/dev/null 2>&1; then
  echo "ERROR: uv not found on PATH."
  echo "       Install: https://docs.astral.sh/uv/getting-started/installation/"
  exit 1
fi

R_AVAILABLE=1
if ! command -v Rscript >/dev/null 2>&1; then
  R_AVAILABLE=0
fi

if command -v quarto >/dev/null 2>&1; then
  echo "Detected: quarto ($(quarto --version 2>/dev/null | head -1))"
else
  echo "Detected: quarto (not installed)"
fi
if [ "$R_AVAILABLE" -eq 1 ]; then
  echo "Detected: Rscript ($(Rscript --version 2>&1 | head -1))"
else
  echo "Detected: Rscript (not installed)"
fi
echo "Detected: uv ($(uv --version 2>&1 | head -1))"
echo

# ── sync python deps ──────────────────────────────────────────────────────
echo "Syncing Python dependencies with uv..."
uv sync
echo

# ── run steps ─────────────────────────────────────────────────────────────
FAILED=0
PASSED=0
STEP=0
TOTAL=6

RED=$'\033[1;31m'
RESET=$'\033[0m'

run_step() {
  STEP=$((STEP + 1))
  local name="$1"; shift
  echo "[$STEP/$TOTAL] $name"
  if "$@"; then
    echo "  PASS: $name"
    PASSED=$((PASSED + 1))
  else
    printf '  %sFAIL: %s%s\n' "$RED" "$name" "$RESET"
    FAILED=$((FAILED + 1))
  fi
  echo
}

run_step "01 Cohort"             uv run python code/01_cohort.py
run_step "02 Wide dataset"       uv run python code/02_wide_dataset.py
run_step "02b Daily SOFA"        uv run python code/02b_sofa.py
run_step "03 SAT"                uv run python code/03_sat.py
run_step "04 SBT"                uv run python code/04_sbt_both.py
run_step "05 Analysis dataset"   uv run python code/05_analysis_dataset.py

if [ "$R_AVAILABLE" -eq 1 ]; then
  run_step "06 R analysis (ABTRISE_run_all.R)" \
    Rscript --vanilla code/ABTRISE_run_all.R
else
  STEP=$((STEP + 1))
  echo "[$STEP/$TOTAL] R analysis"
  echo "  SKIPPED: Rscript not found on PATH."
  echo
  echo "  Run the R analysis from RStudio instead:"
  echo
  echo "    1. Open RStudio."
  echo "    2. Open this file:  code/ABTRISE_run_all.R"
  echo "    3. Click Source (or press Cmd/Ctrl + Shift + Enter)."
  echo
  echo "  That single source() will run all four analysis scripts in order:"
  echo "    - code/ABTRISE_01_setup_c.R          (auto-sourced by the next three)"
  echo "    - code/ABTRISE_02_criterion_c.R"
  echo "    - code/ABTRISE_345_outcomes_c.R"
  echo "    - code/ABTRISE_06_benchmarking_c.R"
  echo
  echo "  First-time R setup (paste into the RStudio console):"
  echo "    install.packages(c(\"here\",\"jsonlite\",\"arrow\",\"dplyr\",\"tidyr\",\"stringr\","
  echo "                       \"ggplot2\",\"patchwork\",\"readr\",\"lme4\",\"glmmTMB\","
  echo "                       \"survival\",\"tidycmprsk\",\"epiR\",\"blandr\",\"splines\","
  echo "                       \"broom\",\"broom.mixed\",\"purrr\",\"scales\",\"forcats\"))"
  echo
  FAILED=$((FAILED + 1))
fi

# ── summary ───────────────────────────────────────────────────────────────
echo "────────────────────────────────────────"
echo "  SUMMARY"
echo "────────────────────────────────────────"
echo "  Passed: $PASSED / $TOTAL"
if [ "$FAILED" -gt 0 ]; then
  printf '  %sFailed: %d / %d%s\n' "$RED" "$FAILED" "$TOTAL" "$RESET"
else
  echo "  Failed: $FAILED / $TOTAL"
fi
echo "  Log:    $LOG_FILE"
echo "  Finished: $(date)"

if [ "$FAILED" -gt 0 ]; then
  exit 1
fi
