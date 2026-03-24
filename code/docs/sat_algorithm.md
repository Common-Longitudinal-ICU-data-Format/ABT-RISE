# SAT Algorithm Documentation

Rule-based detection of Spontaneous Awakening Trial (SAT) eligibility and delivery
from electronic health record (EHR) data, using the CLIF (Common Longitudinal ICU
Format) wide dataset.

**Code**: `03_sat.py` (Marimo notebook)
**Input**: `wide_dataset.parquet` (produced by `02_wide_dataset.py`)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Input Data](#2-input-data)
3. [Calendar-Day Grouping](#3-calendar-day-grouping)
4. [Overnight Screening](#4-overnight-screening)
5. [4-Hour Continuous Eligibility Check](#5-4-hour-continuous-eligibility-check)
6. [Eligible Event Flagging](#6-eligible-event-flagging)
7. [SAT Delivery Detection](#7-sat-delivery-detection)
8. [Day-Level Aggregation](#8-day-level-aggregation)
9. [Concordance Analysis](#9-concordance-analysis)
10. [Hospital Summary](#10-hospital-summary)
11. [Output Files](#11-output-files)
12. [Data Flow Diagram](#12-data-flow-diagram)

---

## 1. Overview

The SAT algorithm answers two questions for every ventilator-day in the ICU:

1. **Was this patient eligible for a Spontaneous Awakening Trial?**
   A patient is eligible when they have been on invasive mechanical ventilation (IMV)
   with active sedation and no paralytics for at least 4 continuous hours during the
   overnight screening window.

2. **Was a SAT actually delivered?**
   Delivery is detected when sedation medications drop to zero (or near-zero) for at
   least 30 minutes while the patient remains on IMV.

The algorithm produces a day-level dataset (`sat_day_level.parquet`) where each row
is one ventilator-day with eligibility status, algorithmic delivery flags, and EHR
ground-truth labels. This dataset is then used for concordance analysis against
nurse-charted SAT flowsheet data.

---

## 2. Input Data

### 2.1 The Wide Dataset (`02_wide_dataset.py`)

The wide dataset merges four data sources onto a common time spine — the union of all
unique `(hospitalization_id, recorded_dttm)` pairs:

| Source | Key columns |
|--------|-------------|
| **Respiratory support** | `device_category`, `mode_category`, `fio2_set`, `peep_set`, etc. |
| **SpO2 vitals** | `spo2` |
| **Medications** | 18 individual drug columns + 4 consolidated dose columns |
| **Assessments** | 10 flowsheet columns (SAT/SBT screens & deliveries, RASS, GCS) |

### 2.2 Key Columns Used by the SAT Algorithm

**Ventilation status:**
- `device_category` — respiratory device; the algorithm checks for `"IMV"` (invasive mechanical ventilation)

**6 sedation medications (individual drug doses):**
- `fentanyl`, `propofol`, `lorazepam`, `midazolam`, `hydromorphone`, `morphine`

**3 paralytic medications:**
- `cisatracurium`, `vecuronium`, `rocuronium`

**4 consolidated dose summary columns (computed in `02_wide_dataset.py`):**

| Column | Definition | Used for |
|--------|-----------|----------|
| `min_sedation_dose` | `min()` across all 6 sedation drugs (nulls skipped) | Stoppage detection — a value of 0 means at least one drug was recorded at dose 0 |
| `min_sedation_dose_2` | `min()` across all 6 sedation drugs **where dose > 0** (zeros treated as null) | Eligibility check — `> 0` means at least one drug is actively being given at a positive dose |
| `min_sedation_dose_non_ops` | `min()` across non-opioid sedatives only (`propofol`, `lorazepam`, `midazolam`), null filled to 0 | Modified delivery detection |
| `max_paralytics` | `max()` across all 3 paralytics, null filled to 0 | Paralytic exclusion — `> 0` means at least one paralytic is active |

**4 EHR ground-truth flowsheet columns:**
- `sat_screen_pass_fail` — nurse-charted SAT screen result
- `sat_screen_performed` — nurse-charted SAT screen performed
- `sat_delivery_pass_fail` — nurse-charted SAT delivery result
- `sat_delivery_performed` — nurse-charted SAT delivery performed

### 2.3 Forward-Fill Policy

After joining sources onto the spine, forward-fill is applied within each patient
(`over("hospitalization_id")`). Crucially, the SAT/SBT flowsheet columns are
**excluded from forward-fill**:

> Forward-filling point-in-time flowsheet flags would bleed a single nurse charting
> event across all subsequent timestamps in the day, falsely inflating the ground
> truth. These columns are left as-is so they only appear at the exact time they were
> charted.

---

## 3. Calendar-Day Grouping

Observations are grouped into ventilator-days using **midnight-to-midnight calendar
days** (`VENT_DAY_ANCHOR_HOUR = 0`).

```
icu_day_date = recorded_dttm truncated to date (midnight)
```

Each patient's vent-days are numbered sequentially:

```
icu_day = dense_rank(icu_day_date) over hospitalization_id
```

A composite key ties each row to its vent-day:

```
hosp_id_icu_day = "{hospitalization_id}_day_{icu_day}"
```

> **Why calendar days?** Calendar-day grouping aligns with clinical workflows (shift
> changes typically happen at fixed clock times) and makes the overnight screening
> window straightforward to define. The tradeoff is that some observations at 22:00-
> 23:59 fall in the *previous* calendar day's group, requiring a cross-day IMV join
> later (see [Section 8](#8-day-level-aggregation)).

---

## 4. Overnight Screening

Before running the expensive 4-hour eligibility loop, each vent-day is quickly
screened using observations in the **22:00-06:00 overnight window**.

### 4.1 Window Definition

For a vent-day with `icu_day_date = D`:
- **Start**: `D - 2 hours` = previous day 22:00
- **End**: `D + 6 hours` = current day 06:00

> **Why 22:00-06:00?** This 8-hour window is designed to capture the overnight period
> when SAT eligibility is clinically assessed. It provides enough time for a 4-hour
> continuous segment to occur while focusing on the pre-rounding period.

### 4.2 Three Quick-Check Conditions

For each `hosp_id_icu_day`, three boolean flags are computed over the overnight window:

| Flag | Condition | Passing value |
|------|-----------|---------------|
| `has_imv` | Any observation with `device_category == "IMV"` | `True` |
| `has_sedation` | Any observation with `min_sedation_dose_2 > 0` | `True` |
| `has_paralytics` | Any observation with `max_paralytics > 0` | `False` (must be absent) |

A day passes screening when: `has_imv AND has_sedation AND NOT has_paralytics`.

### 4.3 Hierarchical Failure Reasons

Days that fail screening receive exactly one reason, assigned in priority order:

1. `"No overnight data"` — no observations in the 22:00-06:00 window at all
2. `"Paralytics present"` — any paralytic was charted during the window
3. `"No IMV"` — patient was not on invasive mechanical ventilation
4. `"No active sedation"` — no sedation drug with a positive dose
5. `"< 4h continuous window"` — passed screening but failed the 4-hour check (default)

### 4.4 CONSORT Counts

The screening produces a stepwise attrition funnel (CONSORT diagram):

```
Total vent-days
  └─ Has overnight data
       └─ No paralytics
            └─ Has IMV
                 └─ Has active sedation
                      └─ >= 4h continuous (computed next)
```

---

## 5. 4-Hour Continuous Eligibility Check

This is the core eligibility algorithm. It runs on every patient-day that passed the
overnight screening.

### 5.1 Row-Level Condition Check

Each row is tagged with a binary flag:

```
all_conditions_check = (device_category == "IMV")
                     & (min_sedation_dose_2 > 0)
                     & (max_paralytics <= 0)
```

This is 1 when the patient is simultaneously on IMV, receiving active sedation, and
not receiving paralytics at that exact timestamp.

### 5.2 Contiguous Segment Identification

Within the 22:00-06:00 window for each vent-day, the algorithm:

1. Extracts the time-ordered observations and their `all_conditions_check` values
2. Detects transitions (where `all_conditions_check` changes from 0 to 1 or vice versa)
3. Assigns a group number to each contiguous run of identical values
4. Filters to groups where `all_conditions_check == 1`

### 5.3 4-Hour Threshold Measurement

For each contiguous "all conditions met" segment:

1. Compute cumulative duration from the first observation in the segment
2. If the segment's total duration (last timestamp minus first) >= 4 hours, the day
   is **eligible**
3. Record `event_time_at_4_hours` — the timestamp of the observation at (or just after)
   the 4-hour mark, found via `np.searchsorted`

> **Why only the first qualifying segment?** Once one 4-hour segment is found, the
> loop breaks (`break`). A patient only needs to be eligible once per day. Taking the
> first segment ensures the eligible event time is as early as possible, maximizing
> the observation window for delivery detection.

---

## 6. Eligible Event Flagging

After the 4-hour check identifies eligible days, two flags are assigned to rows:

### 6.1 `eligible_event`

Set to `1.0` on exactly **one row per eligible day**: the first observation at or
after `event_time_at_4_hours`. This marks the moment eligibility was established and
SAT delivery assessment can begin.

### 6.2 Last-Observation Guard

If the `eligible_event` row happens to be the very last observation for that patient
(across all days), the flag is removed (set to null).

> **Why?** SAT delivery detection requires looking at a 30-minute forward window.
> If there is no subsequent data after the eligible event, we cannot determine whether
> sedation was actually stopped. Flagging such a row as eligible would produce a
> guaranteed false negative.

### 6.3 `on_vent_and_sedation`

Set to `1` on **all rows** belonging to days that have an `eligible_event`. This is
a convenience flag used to filter down to eligible-day rows for the delivery detection
loop.

---

## 7. SAT Delivery Detection

Delivery detection runs only on rows from eligible days (`on_vent_and_sedation == 1`).

### 7.1 `rank_sedation` — Sedation Stoppage Detection

Two rank columns are computed as cumulative counts of zero-dose observations per day:

| Column | Zero-dose column | Drugs included |
|--------|-----------------|----------------|
| `rank_sedation` | `min_sedation_dose == 0` | All 6 (sedatives + opioids) |
| `rank_sedation_non_ops` | `min_sedation_dose_non_ops == 0` | 3 non-opioid sedatives only |

When `rank_sedation` is **not null**, it means `min_sedation_dose == 0` at that row
— i.e., at least one sedation drug was explicitly recorded at dose zero. This is
treated as a "sedation stoppage" event and triggers delivery evaluation.

> **Why `min_sedation_dose` (not `_2`) for stoppage?** `min_sedation_dose` includes
> zero values, so `== 0` captures the moment a drug dose is charted as zero — an
> explicit stoppage signal. `min_sedation_dose_2` excludes zeros by design (treating
> them as null), so it cannot detect stoppages. The two columns serve complementary
> purposes:
>
> - `min_sedation_dose_2 > 0` for eligibility = "is the patient actively sedated?"
> - `min_sedation_dose == 0` for delivery = "did a sedation drug get stopped?"

### 7.2 The 30-Minute Forward Window

At each sedation-stoppage row (non-null `rank_sedation`), the algorithm looks forward
30 minutes and checks:

1. **IMV maintained**: `device_category == "IMV"` for all observations in the forward
   window. If the patient loses ventilation, the row is skipped.
2. **Medication levels**: whether specific drug doses remain at zero for the full
   30-minute window.

### 7.3 Two Delivery Flags

| Flag | Condition | Medications checked |
|------|-----------|-------------------|
| `SAT_primary_delivery` | All 6 sedation meds are zero or null for next 30 min | `fentanyl`, `propofol`, `lorazepam`, `midazolam`, `hydromorphone`, `morphine` |
| `SAT_modified_delivery` | All 3 non-opioid sedatives are zero or null for next 30 min | `propofol`, `lorazepam`, `midazolam` |

- **Primary** = full sedation vacation (all 6 meds off)
- **Modified** = sedative vacation only (opioids may continue for pain management)

A value of `1.0` means delivery was detected at that row; null means it was not.

---

## 8. Day-Level Aggregation

The row-level flags are collapsed into one row per ventilator-day.

### 8.1 Vent-Day Filter

Not all `hosp_id_icu_day` values represent true ventilator-days. Two conditions must
both be met:

1. **Post-intubation**: the vent-day date is on or after the patient's `intubation_time`
   (truncated to date)
2. **IMV in overnight window**: at least one IMV observation exists in the 22:00-06:00
   window for that day

> **Why the cross-day IMV join?** With calendar-day grouping, an observation at 22:30
> belongs to the *previous* calendar day's `hosp_id_icu_day`. But that 22:30 IMV
> observation is clinically relevant to *today's* overnight window (22:00-06:00). The
> cross-day join re-checks IMV presence using the actual clock-time window rather than
> the grouping key, ensuring days aren't incorrectly excluded.

### 8.2 Key Aggregated Columns

| Column | Logic |
|--------|-------|
| `is_eligible` | `True` if `eligible_event == 1` on any row that day |
| `SAT_primary_delivery` | `max()` of the flag across all rows in the day (1 if any row flagged) |
| `SAT_modified_delivery` | Same as above for modified delivery |
| `sat_ground_truth` | `1.0` if **any** of the 4 EHR flowsheet columns is 1 **and** the day was eligible; null otherwise |

### 8.3 Delivery Failure Reasons

For eligible days where delivery was not detected:

| Column | Meaning |
|--------|---------|
| `SAT_primary_delivery_failure` | `"No sedation stoppage"` if `rank_sedation` was never non-null; otherwise `"Meds not zero for 30min or lost IMV"` |
| `SAT_modified_delivery_failure` | `"No non-opioid sedation stoppage"` if `rank_sedation_non_ops` was never non-null; otherwise `"Non-opioid meds not zero for 30min or lost IMV"` |

---

## 9. Concordance Analysis

The algorithmic flags are compared against the EHR ground truth using standard
classification metrics.

### 9.1 Metrics Computed

For each of `SAT_primary_delivery` and `SAT_modified_delivery` vs `sat_ground_truth`:

- **Confusion matrix** (TP, FP, FN, TN)
- **Accuracy**, **Precision**, **Recall**, **F1**, **Specificity**
- **Cohen's Kappa** with 95% bootstrap confidence interval (2000 resamples, seed=42)

### 9.2 Landis-Koch Interpretation Scale

| Kappa range | Interpretation |
|-------------|---------------|
| < 0.00 | Poor |
| 0.00 - 0.20 | Slight |
| 0.21 - 0.40 | Fair |
| 0.41 - 0.60 | Moderate |
| 0.61 - 0.80 | Substantial |
| 0.81 - 1.00 | Almost Perfect |

### 9.3 Outputs

- `delivery_concordance_summary.csv` — all metrics in one table
- `confusion_matrix_{flag_name}.png` — JAMA-style heatmap per flag (cividis colormap)

---

## 10. Hospital Summary

Per-hospital SAT rates for multi-hospital sites (when `hospital_id` is available in
the cohort):

- `eligible_event_count` — number of eligible vent-days
- `{flag}_count` / `pct_{flag}` — count and percentage for each delivery flag and
  ground truth

For single-site cohorts, one summary row is produced for the entire site.

Output: `sat_stats_{SITE_NAME}.csv`

---

## 11. Output Files

All outputs are written to `output_phi/sat_standard/`.

| File | Description |
|------|-------------|
| `sat_day_level.parquet` | Day-level dataset: one row per vent-day with eligibility, delivery flags, ground truth, and failure reasons |
| `consort_sat_{SITE_NAME}.json` | CONSORT attrition funnel as structured JSON |
| `consort_sat_{SITE_NAME}.png` | CONSORT diagram visualization |
| `delivery_concordance_summary.csv` | Concordance metrics (kappa, accuracy, precision, recall, F1, etc.) |
| `confusion_matrix_SAT_primary_delivery.png` | Confusion matrix heatmap for primary delivery |
| `confusion_matrix_SAT_modified_delivery.png` | Confusion matrix heatmap for modified delivery |
| `sat_stats_{SITE_NAME}.csv` | Per-hospital (or site-level) SAT rate summary |

---

## 12. Data Flow Diagram

```
                        ┌───────────────────────────┐
                        │   wide_dataset.parquet     │
                        │  (from 02_wide_dataset.py) │
                        └─────────────┬─────────────┘
                                      │
                                      v
                        ┌─────────────────────────────┐
                        │  Calendar-Day Grouping       │
                        │  icu_day_date, icu_day │
                        │  hosp_id_icu_day             │
                        └─────────────┬───────────────┘
                                      │
                                      v
                    ┌─────────────────────────────────────┐
                    │   Overnight Screening (22:00-06:00)  │
                    │   Quick-check: IMV? Sedation? No     │
                    │   paralytics? Hierarchical failure    │
                    │   reasons + CONSORT counts            │
                    └──────────────┬──────────────────────┘
                                   │
                        passes screening
                                   │
                                   v
                ┌──────────────────────────────────────────┐
                │  4-Hour Continuous Eligibility Check       │
                │  Contiguous segments where:                │
                │    IMV + sedation (dose>0) + no paralytics │
                │  First segment >= 4h → eligible            │
                │  Records event_time_at_4_hours             │
                └──────────────┬───────────────────────────┘
                               │
                        eligible days
                               │
                               v
                ┌──────────────────────────────────────┐
                │   Eligible Event Flagging              │
                │   eligible_event = first obs >= 4h     │
                │   mark on that row (with last-obs      │
                │   guard), on_vent_and_sedation on all  │
                │   rows of eligible days                │
                └──────────────┬───────────────────────┘
                               │
                               v
                ┌──────────────────────────────────────────┐
                │   SAT Delivery Detection                  │
                │   rank_sedation: sedation stoppage events  │
                │   At each stoppage, look forward 30 min:   │
                │     Primary: all 6 meds zero?              │
                │     Modified: 3 non-opioid meds zero?      │
                └──────────────┬───────────────────────────┘
                               │
                               v
                ┌──────────────────────────────────────────┐
                │   Day-Level Aggregation                   │
                │   Vent-day filter (post-intubation +      │
                │   cross-day IMV join)                      │
                │   is_eligible, delivery flags (max),       │
                │   sat_ground_truth, failure reasons         │
                └──────┬──────────────┬───────────────────┘
                       │              │
              ┌────────┘              └──────────┐
              v                                  v
   ┌────────────────────┐           ┌───────────────────────┐
   │  Concordance        │           │  Hospital Summary     │
   │  Cohen's Kappa +    │           │  Per-hospital SAT     │
   │  bootstrap CI,      │           │  rates & counts       │
   │  confusion matrices │           └───────────┬───────────┘
   └────────┬───────────┘                        │
            │                                    v
            v                        sat_stats_{SITE}.csv
   delivery_concordance_summary.csv
   confusion_matrix_*.png

              ┌─────────────────────────────────┐
              │   Final Outputs                  │
              │   sat_day_level.parquet           │
              │   consort_sat_{SITE}.json / .png  │
              └─────────────────────────────────┘
```

---

## Key Design Decisions

### Why `min_sedation_dose_2` for eligibility but `min_sedation_dose` for stoppage?

- **Eligibility** asks: "Is the patient actively receiving sedation?" →
  `min_sedation_dose_2 > 0` answers this because it only considers positive doses
  (zeros are treated as null, so they don't bring the minimum down).
- **Stoppage detection** asks: "Was a sedation drug explicitly charted at zero?" →
  `min_sedation_dose == 0` captures this because it includes zeros in the minimum.

### Why exclude SAT/SBT flowsheet columns from forward-fill?

These are point-in-time charting events. Forward-filling would propagate a single
nurse documentation event across all subsequent timestamps, inflating the ground truth.

### Why is the cross-day IMV join needed?

Calendar-day grouping assigns a 22:30 observation to the previous calendar day's
`hosp_id_icu_day`. But that observation falls within the current day's 22:00-06:00
overnight window. The cross-day join re-checks IMV using clock-time boundaries rather
than grouping keys.

### Why only the first qualifying 4-hour segment per day?

One eligible segment is sufficient to establish SAT eligibility. Taking the earliest
one maximizes the remaining time window available for delivery detection.

### Why the last-observation guard on `eligible_event`?

Delivery detection requires a 30-minute forward window. If the eligible event is the
patient's very last observation, there is no data to evaluate delivery against — the
result would be an automatic false negative. Removing the flag prevents this artifact.
