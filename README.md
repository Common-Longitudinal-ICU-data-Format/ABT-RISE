# ABT-RISE

Awakening & Breathing Trial Rule-based ICU Signature from EHR.

## Objective

ABT-RISE validates a rule-based algorithm that identifies Spontaneous Awakening Trials (SAT) and Spontaneous Breathing Trials (SBT) from routine EHR data, against nurse-charted flowsheet documentation. The pipeline derives an invasive-mechanical-ventilation cohort from CLIF-format ICU data, builds an hourly wide dataset, applies the SAT and SBT detection rules, and runs criterion-, construct-, and benchmarking-level analyses (Lin's CCC, discrete-time logistic, Cox, Fine-Gray competing risks, ZTNB LOS, hospital-level forest plots).

## Required CLIF tables and fields

Version: CLIF 2.x

- `patient` — `patient_id`, `race_category`, `ethnicity_category`, `sex_category`, `death_dttm`
- `hospitalization` — `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `discharge_category`, `age_at_admission`
- `adt` — `hospitalization_id`, `in_dttm`, `out_dttm`, `location_category`
- `respiratory_support` — `hospitalization_id`, `recorded_dttm`, `device_category` (`imv`, `trach collar`), `mode_category`, `tracheostomy`, `fio2_set`, `peep_set`, `resp_rate_set`, `resp_rate_obs`
- `vitals` — `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value`
  - categories used: `spo2`, `weight_kg`, `height_cm`
- `medication_admin_continuous` — `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit`
  - categories used: `norepinephrine`, `epinephrine`, `phenylephrine`, `angiotensin`, `vasopressin`, `dopamine`, `dobutamine`, `milrinone`, `isoproterenol`, `cisatracurium`, `vecuronium`, `rocuronium`, `fentanyl`, `propofol`, `lorazepam`, `midazolam`, `hydromorphone`, `morphine`
- `patient_assessments` — `hospitalization_id`, `recorded_dttm`, `assessment_category`, `numerical_value`, `categorical_value`
  - categories used: `sat_screen_pass_fail`, `sat_screen_performed`, `sat_delivery_pass_fail`, `sat_delivery_performed`, `sbt_screen_pass_fail`, `sbt_screen_performed`, `sbt_delivery_pass_fail`, `sbt_delivery_performed`, `rass`, `gcs_total`
- `hospital_diagnosis` — `hospitalization_id`, `diagnosis_code` (used for Charlson/Elixhauser)

## Cohort identification

Patients on invasive mechanical ventilation (IMV) at any point during a hospitalization, identified from a respiratory-support waterfall over `device_category` in `code/01_cohort.py`. Tracheostomy-only episodes (`device_category = trach collar`) are excluded from the start-of-IMV anchor. The pipeline writes a person-period (hourly) file and a hospitalization-level summary that feed all downstream R analyses.

## Configuration

1. Copy `clif_config_template.json` to `clif_config.json` at the repo root.
2. Edit the four required fields:
   - `site_name` — short identifier for your site (used in output filenames).
   - `data_directory` — absolute path to your CLIF parquet tables.
   - `filetype` — `parquet` (default).
   - `timezone` — IANA name, e.g. `US/Central`.
3. `sofa_batch_size` defaults to `10000` hospitalizations. Lower it to reduce peak memory at the cost of additional source-table reads.
4. Leave `abtrise_input_dir` (default `./output_phi/analysis`) and `abtrise_output_dir` (default `./output_to_share`) at their defaults unless your storage layout differs. Python step 5 writes the analysis parquet files to `abtrise_input_dir`; the R steps read from it and write all shareable tables / figures / model summaries to `abtrise_output_dir` (same folder Python uses for shareable CONSORTs and Table 1s).

### Where outputs land

- `output_phi/` — site-internal intermediate parquet files (in `.gitignore`, do not share).
- `output_to_share/` — everything safe to send to the coordinating center:
  - `consort.{json,png}`, `sat_standard/`, `sbt_both_stabilities/` (from Python).
  - `tables/`, `figures/{a2,a3,a4,a5,a6}/`, `models/{a2,a3,a4,a5,a6}/`, `diagnostics/` (from R).
  - All R output filenames are auto-prefixed with `site_name` (e.g. `RUSH_table1.csv`).

## Prerequisites

- **Python** 3.12+ with [`uv`](https://docs.astral.sh/uv/getting-started/installation/) on PATH.
- **R** 4.x with `Rscript` on PATH (optional — see RStudio fallback below).
- CLIF parquet tables at the path in `data_directory`.

R packages: the pipeline **auto-installs** any missing packages from CRAN on
first run (via `ensure_packages()` in `code/ABTRISE_config.R`), so no manual
step is normally required. If a machine is offline, install them once manually:

```r
install.packages(c("here","jsonlite","arrow","dplyr","tidyr","stringr",
                   "ggplot2","patchwork","readr","lme4","glmmTMB",
                   "survival","tidycmprsk","epiR","blandr",
                   "broom","broom.mixed","purrr","scales","forcats","flextable"))
```

## Running the pipeline

**macOS / Linux:**

```bash
chmod +x run_all.sh
./run_all.sh
```

**Windows (PowerShell):**

```powershell
.\run_all.ps1
```

If PowerShell blocks the script, allow it for the current session only:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
```

Both scripts:

- Sync Python deps via `uv sync`.
- Run the 5 Python steps and the R orchestrator sequentially.
- Tee every line of output to `run_log_YYYYMMDD_HHMM.log` at the repo root.
- Print a final `Passed / Failed` summary and exit non-zero if any step failed.
- **Skip R gracefully** if `Rscript` is not on PATH, and print RStudio fallback instructions.

### RStudio fallback (no Rscript on PATH)

If `Rscript` is missing, the Python steps still run; the wrapper prints these instructions for the R portion:

1. Open RStudio.
2. Open `code/ABTRISE_run_all.R`.
3. Click **Source** (or press Cmd/Ctrl + Shift + Enter).

That single source() runs all four analysis scripts in order:

- `code/ABTRISE_00_setup_c.R` — auto-sourced by the next three
- `code/ABTRISE_02_criterion_c.R`
- `code/ABTRISE_345_outcomes_c.R`
- `code/ABTRISE_06_benchmarking_c.R`

## Pipeline steps

| # | Script | Language | Description |
|---|---|---|---|
| 1 | `code/01_cohort.py` | Python | Load CLIF tables via clifpy, build IMV cohort, write `output_phi/cohort.parquet` + helper tables. |
| 2 | `code/02_wide_dataset.py` | Python | Join cohort + vitals + meds + assessments into hourly wide dataset. |
| 3 | `code/03_sat.py` | Python | Apply SAT rule (standard variant), produce day-level SAT table + Table 1. |
| 4 | `code/04_sbt_both.py` | Python | Apply SBT rule (both stability variants), produce day-level SBT table + Table 1. |
| 5 | `code/05_analysis_dataset.py` | Python | Assemble `file1_person_period.parquet` and `file2_hospitalization_level.parquet`. |
| 6 | `code/ABTRISE_run_all.R` | R | Orchestrator: runs setup, criterion validity (A2), outcomes (A3–A5, SA), and hospital benchmarking (A6). |

## Project structure

```
.
├── code/
│   ├── 01_cohort.py                  # Step 1: IMV cohort
│   ├── 02_wide_dataset.py            # Step 2: hourly wide dataset
│   ├── 03_sat.py                     # Step 3: SAT detection
│   ├── 04_sbt_both.py                # Step 4: SBT detection
│   ├── 05_analysis_dataset.py        # Step 5: person-period + hosp-level files
│   ├── ABTRISE_config.R              # R config loader (reads clif_config.json)
│   ├── ABTRISE_00_setup_c.R          # R: data load + Table 1 (auto-sourced)
│   ├── ABTRISE_02_criterion_c.R     # R: criterion validity (A2)
│   ├── ABTRISE_345_outcomes_c.R     # R: construct validity (A3, A4, A5, SA)
│   ├── ABTRISE_06_benchmarking_c.R  # R: hospital benchmarking (A6)
│   ├── ABTRISE_run_all.R             # R orchestrator
│   ├── helper.py                     # shared Python utilities
│   └── extra/                        # auxiliary notebooks
├── clif_config_template.json         # config template (copy to clif_config.json)
├── run_all.sh                        # pipeline runner (macOS/Linux)
├── run_all.ps1                       # pipeline runner (Windows)
├── pyproject.toml                    # Python deps (uv)
└── README.md
```
