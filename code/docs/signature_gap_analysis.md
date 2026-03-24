# Gap Analysis: ABT-RISE vs CLIF Signature Repo

Comparison of the current ABT-RISE site implementation (`code/`) against the
multi-site signature repo (`CLIF_rule_based_SAT_SBT_signature/`) requirements
defined in `utils/definitions_source_of_truth.py` and the Statistical Analysis
Plan (SAP).

---

## 1. Vent-Day Anchor

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Anchor hour** | `VENT_DAY_ANCHOR_HOUR = 0` (midnight) | `VENT_DAY_ANCHOR_HOUR = 6` (06:00) |
| **File** | `02_wide_dataset.py:330` | `definitions_source_of_truth.py:35` |

**Impact:** Every vent-day boundary, overnight window calculation, and
day-level aggregation shifts by 6 hours. A "day" becomes 06:00 → 05:59 next
day instead of 00:00 → 23:59.

**What needs to change:**
- `02_wide_dataset.py`: Change `VENT_DAY_ANCHOR_HOUR` from `0` to `6`.
- All downstream scripts (`03_sat.py`, `04_sbt*.py`) inherit `icu_day_date`
  from the wide dataset, so no changes needed there for anchoring itself.
- Verify the overnight window math still holds: the 22:00–06:00 window
  relative to a 06:00-anchored day means 22:00 on the *prior calendar day*
  to 06:00 on the *index calendar day* — which is the same clock range but
  now sits differently relative to the vent-day boundary.

---

## 2. SAT RASS Agitation Check

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **RASS check** | Not implemented | RASS < 2 (exclude if agitated) |
| **Threshold** | — | `RASS_AGITATION_THRESHOLD = 2` |
| **File** | `03_sat.py:153–159` | `definitions_source_of_truth.py:74` |

**Current ABT-RISE `all_conditions_check`:**
```python
(device == "imv") & (sedation > 0) & (paralytics <= 0)
```

**Signature requires:**
```python
(device == "imv") & (sedation > 0) & (paralytics <= 0) & (RASS < 2)
```

**What needs to change:**
- `03_sat.py`: Add `& (pl.col("rass").fill_null(0.0) < 2)` to the
  `all_conditions_check` expression in `process_cohort`.
- Add "RASS >= 2 (agitated)" to the overnight screening failure reasons.
- Add `has_rass_agitated` to `failure_flags_df` for UpSet plots.
- Decide on null-handling strategy for RASS: signature fills null → 0
  (assumes non-agitated if not charted).

---

## 3. SBT Controlled Mode Duration

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Required duration** | 6 hours continuous IMV | 12 hours continuous controlled mode |
| **Mode check** | Any IMV | Specific controlled modes only |
| **File** | `04_sbt*.py:167` | `definitions_source_of_truth.py:109` |

**Current ABT-RISE eligibility:**
```python
(device == "imv") & (paralytics <= 0)   # 6h continuous
```

**Signature requires:**
```python
(device == "imv") & (mode in CONTROLLED_MODES) & (paralytics <= 0)   # 12h continuous
```

Where `SBT_CONTROLLED_MODES` = `["assist control-volume control",
"pressure control", "pressure-regulated volume control", "simv"]`

**What needs to change (all 4 SBT files):**
- `04_sbt.py`, `04_sbt_hemo.py`, `04_sbt_resp.py`, `04_sbt_both.py`:
  - Change `_6h_us` → `_12h_us` (12 * 3600 * 1e6).
  - Add `mode_category` check to `all_conditions_check`:
    `& (pl.col("mode_category").str.to_lowercase().is_in(CONTROLLED_MODES))`
  - Update all variable names (`event_time_at_6_hours` → `event_time_at_12_hours`).
  - Update docstrings, comments, CONSORT descriptions.
  - Update failure reason string `"< 6h continuous window"` →
    `"< 12h continuous controlled mode"`.

---

## 4. SBT Hemodynamic Stability

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **NEE threshold** | NEE <= 0.2 | NEE <= 0.2 |
| **Dopamine** | Not checked | Dopamine <= 5.0 mcg/kg/min |
| **Dobutamine** | Not checked | Dobutamine <= 5.0 mcg/kg/min |
| **Vasopressin** | Included in NEE | Any dose allowed (np.inf) |
| **Milrinone** | Not checked | Any dose allowed (np.inf) |
| **Stability duration** | Row-level (in contiguous segment) | >= 2 consecutive hours |
| **File** | `04_sbt_hemo.py:78–81` | `definitions_source_of_truth.py:117–127` |

**Current ABT-RISE:**
```python
hemodynamic_stable = (nee.fill_null(0) <= 0.2)
```

**Signature requires:**
```python
hemodynamic_stable = (
    (nee.fill_null(0) <= 0.2)
    & (dopamine.fill_null(0) <= 5.0)
    & (dobutamine.fill_null(0) <= 5.0)
)
```

**What needs to change:**
- `04_sbt_hemo.py` and `04_sbt_both.py`: Add dopamine and dobutamine
  threshold checks to the `hemodynamic_stable` flag.
- Consider adding a 2h consecutive stability check within the eligibility
  window (see Gap #3b below).

### 3b. Stability Duration (2h consecutive)

The signature repo (`pySBT.py`) requires respiratory and hemodynamic
stability to be sustained for **>= 2 consecutive hours** within the
22:00–06:00 window. ABT-RISE checks stability row-by-row within the
contiguous eligibility segment but does not enforce a separate 2h minimum.

**What needs to change:**
- After computing the per-row stability flag, add a contiguous-segment
  analysis (similar to the 4h/12h eligibility check) that verifies at least
  one 2h block of consecutive stability exists within the overnight window.

---

## 5. SBT Modified Delivery Duration

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Primary (EHR)** | >= 2 min FLIP | >= 2 min FLIP |
| **Modified** | >= 30 min FLIP | >= 5 min FLIP |
| **Extended** | — | >= 30 min FLIP (sensitivity) |
| **File** | `04_sbt.py:313,358–359` | `definitions_source_of_truth.py:131,144` |

**Current flag columns:** `SBT_delivery_2min`, `SBT_delivery_30min`

**Signature expects:** `EHR_Delivery_2mins`, `EHR_Delivery_5mins`,
`EHR_Delivery_30mins`

**What needs to change (all 4 SBT files):**
- Add a 5-minute FLIP check (`_delta_5min_us = int(5 * 60 * 1e6)`).
- Rename columns to match signature schema:
  - `SBT_delivery_2min` → `EHR_Delivery_2mins`
  - New: `EHR_Delivery_5mins`
  - `SBT_delivery_30min` → `EHR_Delivery_30mins`
- Update `_FLAG_COLS`, delivery failure reason columns, CONSORT printouts,
  UpSet plots, concordance analysis, and hospital summary stats.

---

## 6. IMV Episode Classification

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Episode logic** | Not implemented | 72h gap → new episode |
| **File** | — | `definitions_source_of_truth.py:38,253–288` |

The signature repo classifies IMV observations into distinct episodes
separated by >= 72 hours without ventilator support. Episode IDs are used
for:
- Cluster-bootstrap CIs (clustering at episode level)
- Per-episode outcome calculations
- Minimum 6h IMV per episode for inclusion

**What needs to change:**
- `01_cohort.py`: After intubation/extubation detection, apply
  `classify_imv_episodes()` logic — sort IMV observations by timestamp,
  flag gaps > 72h as new episodes, assign `imv_episode_id`.
- Propagate `imv_episode_id` through the cohort parquet and downstream
  wide dataset.
- Exclude episodes with < 6h total IMV time (`MIN_IMV_HOURS = 6`).

---

## 7. Encounter Stitching

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Stitching** | Not implemented | Admissions within 6h → encounter_block |
| **File** | — | `00_cohort_id.py` |
| **Output** | — | `hospitalization_to_block_df.csv` |

The signature repo groups multiple hospitalizations for the same patient
into contiguous "encounter blocks" when discharge-to-readmission gap is
< 6 hours (or configurable). This prevents artificial episode breaks at
administrative boundaries.

**What needs to change:**
- `01_cohort.py`: After loading hospitalizations, sort by
  `(patient_id, admission_dttm)`, compute gap to previous discharge, flag
  gaps < 6h, assign cumulative `encounter_block` IDs.
- Output a mapping file `hospitalization_to_block_df.csv`.
- Use `encounter_block` as the longitudinal unit for episode classification.

---

## 8. Extubation Outcome Definition

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Definition** | 2+ consecutive non-IMV observations | Same |
| **Status** | Already implemented | Matches |
| **File** | `01_cohort.py:393–398` | `definitions_source_of_truth.py:153` |

**No changes needed** — ABT-RISE already implements the 2-row non-IMV
extubation detection at `01_cohort.py:393–398`.

---

## 9. Ventilator-Free Days (VFD-28)

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Computed** | No | Yes (VFD_MAX_DAYS = 28) |
| **File** | — | `definitions_source_of_truth.py:151` |

VFD-28 = days alive and free of mechanical ventilation through hospital
day 28. If patient dies before day 28, VFD-28 = 0.

**What needs to change:**
- `01_cohort.py` or a new outcome computation cell: For each
  hospitalization, compute total ventilator days (from intubation to
  extubation), then `VFD-28 = max(0, 28 - total_vent_days)`. If
  `discharge_category == "Expired"` and death occurs before day 28,
  set VFD-28 = 0.
- Add `vfd_28` and `total_vent_days` columns to cohort output.

---

## 10. Comorbidity Index

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Index** | Charlson (CCI) via `clifpy.calculate_cci` | Elixhauser van Walraven |
| **File** | `03_sat.py:618–626`, `04_sbt.py:769–777` | `utils/elixhauser.py` |

The SAP specifies Elixhauser van Walraven comorbidity score for covariate
adjustment in outcome models.

**What needs to change:**
- Both `03_sat.py` and `04_sbt.py` (Table 1 / covariate sections):
  Replace `calculate_cci` with an Elixhauser implementation.
- The signature repo provides `utils/elixhauser.py` — either import it
  or use an equivalent (e.g., `comorbidipy` package or manual ICD-10
  mapping).
- Rename `cci_score` → `elixhauser_score` (or `elix_vanwalraven`) in
  Table 1 output.
- **Option:** Keep CCI as a secondary metric and add Elixhauser alongside.

---

## 11. Cluster-Bootstrap CIs

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Method** | Simple row-level bootstrap (2000 reps) | BCa cluster-bootstrap at episode level (1000 reps) |
| **Cluster unit** | — | IMV episode |
| **File** | `03_sat.py:986–998`, `04_sbt.py:660–672` | SAP Section 2.7 |

The SAP requires bias-corrected and accelerated (BCa) bootstrap confidence
intervals, with resampling at the **episode level** (not row level) to
account for within-episode correlation.

**What needs to change:**
- `03_sat.py` and `04_sbt.py` (`_kappa_ci_bootstrap` function):
  - Change from row-level to episode-level resampling: resample
    `imv_episode_id`s with replacement, then include all vent-days within
    each sampled episode.
  - Switch from percentile CI to BCa CI (requires computing acceleration
    and bias-correction factors).
  - Change `n_boot` from 2000 to 1000.
- Requires `imv_episode_id` to be available in the day-level dataset
  (depends on Gap #6).

---

## 12. Output Schema Alignment

### Column Renaming

| ABT-RISE Column | Signature Column | File |
|---|---|---|
| `SAT_primary_delivery` | `SAT_EHR_delivery` | `03_sat.py` |
| `SAT_modified_delivery` | `SAT_modified_delivery` | (matches) |
| `SBT_delivery_2min` | `EHR_Delivery_2mins` | `04_sbt*.py` |
| (missing) | `EHR_Delivery_5mins` | `04_sbt*.py` |
| `SBT_delivery_30min` | `EHR_Delivery_30mins` | `04_sbt*.py` |
| `on_vent_and_sedation` | `on_vent_and_sedation` | (matches) |
| `on_vent_eligible` | `eligible_day` | `04_sbt*.py` |
| (missing) | `extubated` | `04_sbt*.py` |
| (missing) | `died` | `04_sbt*.py` |
| (missing) | `vfd_28` | `04_sbt*.py` |
| (missing) | `total_vent_days` | `04_sbt*.py` |
| (missing) | `encounter_block` | All |
| (missing) | `hospital_id` | SAT/SBT outputs |
| (missing) | `patient_id` | SAT/SBT outputs |

### Missing Output Files

| File | Description | Status |
|---|---|---|
| `study_cohort.parquet` | Row-level cohort with episodes | Partially exists as `cohort.parquet` — needs episode IDs, encounter blocks |
| `final_df_SAT.csv` | Day-level SAT with all covariates | Exists as `sat_day_level.parquet` — needs column renames + missing columns |
| `final_df_SBT.csv` | Day-level SBT with outcomes | Exists as `sbt_day_level.parquet` — needs outcomes + column renames |
| `hospitalization_to_block_df.csv` | Encounter stitching map | Does not exist |

---

## 13. Study Period Filter

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Period** | 2018–2024 | 2022–2024 |
| **File** | `01_cohort.py:89–91` | `definitions_source_of_truth.py:27–28` |

**What needs to change:**
- `01_cohort.py`: Narrow admission filter from `2018–2024` to `2022–2024`.

---

## 14. Age Cap

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Min age** | 18 | 18 |
| **Max age** | Not enforced | 119 (CLIF convention) |
| **File** | `01_cohort.py:75` | `definitions_source_of_truth.py:30` |

**What needs to change:**
- `01_cohort.py`: Add `& (pl.col("age_at_admission") <= 119)` to the
  adult filter (line 75).

---

## 15. Tracheostomy Exclusion

| | ABT-RISE | Signature (SAP) |
|---|---|---|
| **Method** | First 24h of admission + trach collar | Pre- and post-waterfall exclusion |
| **File** | `01_cohort.py:236–268` | `00_cohort_id.py` |

The signature repo applies tracheostomy exclusion both before and after
the respiratory support waterfall fill, to catch cases where waterfall
propagation fills in trach status.

**What needs to change:**
- `01_cohort.py`: Add a second tracheostomy check after waterfall
  processing (`resp_wf`) to catch trach cases introduced by forward-fill.

---

## Priority Order for Implementation

| Priority | Gap | Effort | Impact |
|---|---|---|---|
| 1 | Vent-day anchor (06:00) | Low | High — changes all day boundaries |
| 2 | SBT 12h controlled mode | Medium | High — fundamental eligibility change |
| 3 | SBT 5min delivery tier | Low | Medium — adds required output column |
| 4 | SAT RASS check | Low | Medium — adds eligibility criterion |
| 5 | SBT hemodynamic stability | Low | Medium — adds dopamine/dobutamine |
| 6 | Study period 2022–2024 | Low | Medium — narrows cohort |
| 7 | IMV episode classification | Medium | High — needed for bootstrap CIs |
| 8 | VFD-28 outcome | Low | Medium — needed for outcome models |
| 9 | Encounter stitching | Medium | Medium — needed for output schema |
| 10 | Elixhauser comorbidity | Medium | Low — covariate swap |
| 11 | Output column renaming | Low | Low — cosmetic alignment |
| 12 | Cluster-bootstrap CIs | High | Medium — statistical method change |
| 13 | Age cap (119) | Low | Low — likely no effect |
| 14 | Post-waterfall trach exclusion | Low | Low — edge case |
| 15 | 2h stability duration check | Medium | Medium — adds strictness |
