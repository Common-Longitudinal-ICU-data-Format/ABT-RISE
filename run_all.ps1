# ABT-RISE pipeline runner (Windows PowerShell)
#
# Runs 5 Python steps then the R orchestrator. Tees everything to
# run_log_YYYYMMDD_HHMM.log at the repo root. Skips R gracefully and
# prints RStudio fallback instructions if Rscript is not on PATH.
#
# Usage:  .\run_all.ps1
# If blocked by execution policy, run once in the same shell:
#   Set-ExecutionPolicy -Scope Process Bypass

$ProjectRoot = $PSScriptRoot
Set-Location $ProjectRoot

$Timestamp = Get-Date -Format "yyyyMMdd_HHmm"
$LogFile = Join-Path $ProjectRoot "run_log_$Timestamp.log"

function Write-Log {
    param([string]$Message = "")
    $Message | Tee-Object -FilePath $LogFile -Append
}

$RED   = [char]27 + "[1;31m"
$RESET = [char]27 + "[0m"

function Invoke-Step {
    param(
        [string]$Name,
        [scriptblock]$Action
    )
    $script:Step++
    Write-Log "[$($script:Step)/$($script:Total)] $Name"
    try {
        & $Action 2>&1 | Tee-Object -FilePath $LogFile -Append
        if ($LASTEXITCODE -ne $null -and $LASTEXITCODE -ne 0) {
            throw "Exit code $LASTEXITCODE"
        }
        Write-Log "  PASS: $Name"
        $script:Passed++
    } catch {
        Write-Log "  ${RED}FAIL: $Name ($_)${RESET}"
        $script:Failed++
    }
    Write-Log ""
}

Write-Log "ABT-RISE pipeline"
Write-Log "Started: $(Get-Date)"
Write-Log "Log:     $LogFile"
Write-Log ""

# ── config check ──────────────────────────────────────────────────────────
if (-not (Test-Path (Join-Path $ProjectRoot "clif_config.json"))) {
    Write-Log "ERROR: clif_config.json not found at repo root."
    Write-Log "       Copy clif_config_template.json to clif_config.json and edit it."
    exit 1
}

# ── tool checks ───────────────────────────────────────────────────────────
$uvCmd = Get-Command uv -ErrorAction SilentlyContinue
if (-not $uvCmd) {
    Write-Log "ERROR: uv not found on PATH."
    Write-Log "       Install: https://docs.astral.sh/uv/getting-started/installation/"
    exit 1
}

$rscriptCmd = Get-Command Rscript -ErrorAction SilentlyContinue
$RAvailable = [bool]$rscriptCmd

$quartoCmd = Get-Command quarto -ErrorAction SilentlyContinue
if ($quartoCmd) {
    $qv = (& quarto --version 2>$null | Select-Object -First 1)
    Write-Log "Detected: quarto ($qv)"
} else {
    Write-Log "Detected: quarto (not installed)"
}
if ($RAvailable) {
    $rv = (& Rscript --version 2>&1 | Select-Object -First 1)
    Write-Log "Detected: Rscript ($rv)"
} else {
    Write-Log "Detected: Rscript (not installed)"
}
$uvv = (& uv --version 2>&1 | Select-Object -First 1)
Write-Log "Detected: uv ($uvv)"
Write-Log ""

# ── sync python deps ──────────────────────────────────────────────────────
Write-Log "Syncing Python dependencies with uv..."
& uv sync 2>&1 | Tee-Object -FilePath $LogFile -Append
Write-Log ""

# ── run steps ─────────────────────────────────────────────────────────────
$script:Step = 0
$script:Passed = 0
$script:Failed = 0
$script:Total = 6

Invoke-Step "01 Cohort"           { uv run python code/01_cohort.py }
Invoke-Step "02 Wide dataset"     { uv run python code/02_wide_dataset.py }
Invoke-Step "03 SAT"              { uv run python code/03_sat.py }
Invoke-Step "04 SBT"              { uv run python code/04_sbt_both.py }
Invoke-Step "05 Analysis dataset" { uv run python code/05_analysis_dataset.py }

if ($RAvailable) {
    Invoke-Step "06 R analysis (ABTRISE_run_all.R)" {
        Rscript --vanilla code/ABTRISE_run_all.R
    }
} else {
    $script:Step++
    Write-Log "[$($script:Step)/$($script:Total)] R analysis"
    Write-Log "  SKIPPED: Rscript not found on PATH."
    Write-Log ""
    Write-Log "  Run the R analysis from RStudio instead:"
    Write-Log ""
    Write-Log "    1. Open RStudio."
    Write-Log "    2. Open this file:  code/ABTRISE_run_all.R"
    Write-Log "    3. Click Source (or press Ctrl + Shift + Enter)."
    Write-Log ""
    Write-Log "  That single source() will run all four analysis scripts in order:"
    Write-Log "    - code/ABTRISE_01_setup_c.R          (auto-sourced by the next three)"
    Write-Log "    - code/ABTRISE_02_criterion_c.R"
    Write-Log "    - code/ABTRISE_345_outcomes_c.R"
    Write-Log "    - code/ABTRISE_06_benchmarking_c.R"
    Write-Log ""
    Write-Log "  First-time R setup (paste into the RStudio console):"
    Write-Log "    install.packages(c(`"here`",`"jsonlite`",`"arrow`",`"dplyr`",`"tidyr`",`"stringr`","
    Write-Log "                       `"ggplot2`",`"patchwork`",`"readr`",`"lme4`",`"glmmTMB`","
    Write-Log "                       `"survival`",`"tidycmprsk`",`"epiR`",`"blandr`",`"splines`","
    Write-Log "                       `"broom`",`"broom.mixed`",`"purrr`",`"scales`",`"forcats`"))"
    Write-Log ""
    $script:Failed++
}

# ── summary ───────────────────────────────────────────────────────────────
Write-Log "----------------------------------------"
Write-Log "  SUMMARY"
Write-Log "----------------------------------------"
Write-Log "  Passed: $($script:Passed) / $($script:Total)"
if ($script:Failed -gt 0) {
    Write-Log "  ${RED}Failed: $($script:Failed) / $($script:Total)${RESET}"
} else {
    Write-Log "  Failed: $($script:Failed) / $($script:Total)"
}
Write-Log "  Log:    $LogFile"
Write-Log "  Finished: $(Get-Date)"

if ($script:Failed -gt 0) {
    exit 1
}
