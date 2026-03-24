# SBT Algorithm Documentation

Rule-based detection of Spontaneous Breathing Trial (SBT) eligibility and delivery
from electronic health record (EHR) data, using the CLIF (Common Longitudinal ICU
Format) wide dataset.

**Code**: `04_sbt.py` (standard), `04_sbt_hemo.py` (hemodynamic), `04_sbt_resp.py` (respiratory), `04_sbt_both.py` (combined) — all Marimo notebooks
**Input**: `wide_dataset.parquet` (produced by `02_wide_dataset.py`)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Input Data](#2-input-data)
3. [Calendar-Day Grouping](#3-calendar-day-grouping)
4. [Overnight Screening (22:00–06:00)](#4-overnight-screening-220006-00)
5. [6-Hour Continuous Eligibility Check](#5-6-hour-continuous-eligibility-check)
6. [Stability Variants](#6-stability-variants)
7. [Eligible Event Flagging](#7-eligible-event-flagging)
8. [SBT Delivery Detection](#8-sbt-delivery-detection)
9. [Day-Level Aggregation](#9-day-level-aggregation)
10. [Concordance Analysis](#10-concordance-analysis)
11. [UpSet / Failure Analysis](#11-upset--failure-analysis)
12. [Output Files](#12-output-files)
13. [Data Flow Diagram](#13-data-flow-diagram)
14. [Key Design Decisions](#14-key-design-decisions)

---

## 1. Overview

The SBT algorithm answers two questions for every ventilator-day in the ICU:

1. **Was this patient eligible for a Spontaneous Breathing Trial?**
   A patient is eligible when they have been on invasive mechanical ventilation (IMV)
   with no paralytics for at least 6 continuous hours during the overnight screening
   window. Unlike SAT, SBT does **not** require active sedation.

2. **Was an SBT actually delivered?**
   Delivery is detected when the ventilator mode transitions ("FLIPs") to pressure
   support, CPAP, or T-piece with low support settings (`pressure_support_set <= 8` and
   `peep_set <= 8`) and this mode is sustained for at least 2 or 30 minutes.

The algorithm produces a day-level dataset (`sbt_day_level.parquet`) where each row
is one ventilator-day with eligibility status, algorithmic delivery flags, and EHR
ground-truth labels. This dataset is then used for concordance analysis against
nurse-charted SBT flowsheet data.

Four variants of the algorithm exist, differing only in which stability conditions
are added to the 6-hour eligibility check:

| Variant | Code | Output directory | Extra eligibility condition |
|---------|------|-----------------|---------------------------|
| **Standard** | `04_sbt.py` | `sbt_standard/` | None |
| **Hemodynamic** | `04_sbt_hemo.py` | `sbt_hemodynamic/` | NEE ≤ 0.2 |
| **Respiratory** | `04_sbt_resp.py` | `sbt_respiratory/` | FiO2 ≤ 0.5 & PEEP ≤ 8 & SpO2 ≥ 88 |
| **Both** | `04_sbt_both.py` | `sbt_both_stabilities/` | Hemodynamic + Respiratory |

---

## 2. Input Data

### 2.1 The Wide Dataset (`02_wide_dataset.py`)

The wide dataset merges four data sources onto a common time spine — the union of all
unique `(hospitalization_id, recorded_dttm)` pairs:

| Source | Key columns |
|--------|-------------|
| **Respiratory support** | `device_category`, `mode_category`, `mode_name`, `fio2_set`, `peep_set`, `pressure_support_set` |
| **SpO2 vitals** | `spo2` |
| **Medications** | 18 individual drug columns + 4 consolidated dose columns |
| **Assessments** | 10 flowsheet columns (SAT/SBT screens & deliveries, RASS, GCS) |

### 2.2 Key Columns Used by the SBT Algorithm

**Ventilation status:**
- `device_category` — respiratory device; the algorithm checks for `"IMV"` (invasive mechanical ventilation)
- `mode_category` — ventilator mode; used in FLIP detection (checks for "pressure support" or "cpap")
- `mode_name` — ventilator mode name; used to detect T-piece (`^t-?piece$`)
- `pressure_support_set` — pressure support level; FLIP requires `<= 8`
- `peep_set` — PEEP level; used in FLIP detection (`<= 8`) and respiratory stability (`<= 8`)

**Paralytic medications:**
- `max_paralytics` — `max()` across `cisatracurium`, `vecuronium`, `rocuronium`, null filled to 0. Value `> 0` means at least one paralytic is active.

**Stability columns (variant-dependent):**
- `nee` — norepinephrine equivalent (hemodynamic stability)
- `fio2_set` — fraction of inspired oxygen (respiratory stability)
- `spo2` — oxygen saturation (respiratory stability)

**4 EHR ground-truth flowsheet columns:**
- `sbt_screen_pass_fail` — nurse-charted SBT screen result
- `sbt_screen_performed` — nurse-charted SBT screen performed
- `sbt_delivery_pass_fail` — nurse-charted SBT delivery result
- `sbt_delivery_performed` — nurse-charted SBT delivery performed

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
> later (see [Section 9](#9-day-level-aggregation)).

---

## 4. Overnight Screening (22:00–06:00)

Before running the expensive 6-hour eligibility loop, each vent-day is quickly
screened using observations in the **22:00-06:00 overnight window**.

### 4.1 Window Definition

For a vent-day with `icu_day_date = D`:
- **Start**: `D - 2 hours` = previous day 22:00
- **End**: `D + 6 hours` = current day 06:00

> **Why 22:00-06:00?** This 8-hour window is designed to capture the overnight period
> when SBT eligibility is clinically assessed. It provides enough time for a 6-hour
> continuous segment to occur while focusing on the pre-rounding period.

### 4.2 Three Quick-Check Conditions

For each `hosp_id_icu_day`, three boolean flags are computed over the overnight window:

| Flag | Condition | Passing value |
|------|-----------|---------------|
| `has_overnight_data` | Any observation exists in the 22:00-06:00 window | `True` |
| `has_imv` | Any observation with `device_category == "IMV"` | `True` |
| `has_paralytics` | Any observation with `max_paralytics > 0` | `False` (must be absent) |

A day passes screening when: `has_overnight_data AND has_imv AND NOT has_paralytics`.

> **Key difference from SAT:** SBT screening does **not** check for active sedation.
> SAT requires `has_sedation` (any `min_sedation_dose_2 > 0`), but SBT does not —
> a patient can be eligible for SBT regardless of sedation status.

### 4.3 Hierarchical Failure Reasons

Days that fail screening receive exactly one reason, assigned in priority order:

1. `"No overnight data"` — no observations in the 22:00-06:00 window at all
2. `"No IMV"` — patient was not on invasive mechanical ventilation
3. `"Paralytics present"` — any paralytic was charted during the window
4. `"< 6h continuous window"` — passed screening but failed the 6-hour check (default)

### 4.4 CONSORT Counts

The screening produces a stepwise attrition funnel (CONSORT diagram):

```
Total vent-days
  └─ Has overnight data
       └─ Has IMV
            └─ No paralytics (implied by >= 6h continuous)
                 └─ >= 6h continuous (computed next)
```

---

## 5. 6-Hour Continuous Eligibility Check

This is the core eligibility algorithm. It runs on every patient-day that passed the
overnight screening.

### 5.1 Row-Level Condition Check (Standard)

Each row is tagged with a binary flag:

```
all_conditions_check = (device_category == "IMV")
                     & (max_paralytics <= 0)
```

This is 1 when the patient is simultaneously on IMV and not receiving paralytics at
that exact timestamp. Unlike SAT's 4-hour check, there is no sedation requirement.

### 5.2 Contiguous Segment Identification

Within the 22:00-06:00 window for each vent-day, the algorithm:

1. Extracts the time-ordered observations and their `all_conditions_check` values
2. Detects transitions (where `all_conditions_check` changes from 0 to 1 or vice versa)
3. Assigns a group number to each contiguous run of identical values
4. Filters to groups where `all_conditions_check == 1`

### 5.3 6-Hour Threshold Measurement

For each contiguous "all conditions met" segment:

1. Compute cumulative duration from the first observation in the segment
2. If the segment's total duration (last timestamp minus first) >= 6 hours, the day
   is **eligible**
3. Record `event_time_at_6_hours` — the timestamp of the observation at (or just after)
   the 6-hour mark, found via `np.searchsorted`

> **Why only the first qualifying segment?** Once one 6-hour segment is found, the
> loop breaks (`break`). A patient only needs to be eligible once per day. Taking the
> first segment ensures the eligible event time is as early as possible, maximizing
> the observation window for delivery detection.

---

## 6. Stability Variants

The four SBT scripts share identical logic for overnight screening, contiguous segment
detection, FLIP delivery detection, and day-level aggregation. They differ only in
what conditions are included in `all_conditions_check` during the 6-hour eligibility
loop, and what additional quick-check flags appear in the overnight screening.

### 6.1 Standard (`04_sbt.py`)

No additional stability conditions. The base eligibility check is:

```
all_conditions_check = (device_category == "IMV")
                     & (max_paralytics <= 0)
```

**Overnight quick-check flags:** `has_paralytics`
**Failure UpSet columns:** `has_paralytics`, `no_6h_window`

### 6.2 Hemodynamic Stability (`04_sbt_hemo.py`)

Adds a hemodynamic stability requirement based on norepinephrine equivalent (NEE):

```python
hemodynamic_stable = (nee.fill_null(0.0) <= 0.2).cast(Int32)
```

The row-level eligibility becomes:

```
all_conditions_check = (device_category == "IMV")
                     & (max_paralytics <= 0)
                     & (hemodynamic_stable == 1)
```

**Null handling:** `nee` nulls are filled with `0.0` (assumes no vasopressor = stable).

**Overnight quick-check flags:** `has_paralytics`, `has_hemo_stable`
**Failure UpSet columns:** `has_paralytics`, `hemo_unstable`, `no_6h_window`

### 6.3 Respiratory Stability (`04_sbt_resp.py`)

Adds a respiratory stability requirement based on FiO2, PEEP, and SpO2:

```python
respiratory_stable = (
    (fio2_set.fill_null(0.21) <= 0.5)
    & (peep_set.fill_null(0.0) <= 8)
    & (spo2.fill_null(100.0) >= 88)
).cast(Int32)
```

The row-level eligibility becomes:

```
all_conditions_check = (device_category == "IMV")
                     & (max_paralytics <= 0)
                     & (respiratory_stable == 1)
```

**Null handling:**
| Column | Null fill value | Rationale |
|--------|----------------|-----------|
| `fio2_set` | `0.21` | Room air (most conservative assumption) |
| `peep_set` | `0.0` | No PEEP (most conservative assumption) |
| `spo2` | `100.0` | Perfect saturation (most conservative assumption) |

**Overnight quick-check flags:** `has_paralytics`, `has_resp_stable`
**Failure UpSet columns:** `has_paralytics`, `resp_unstable`, `no_6h_window`

### 6.4 Both Stabilities (`04_sbt_both.py`)

Requires hemodynamic AND respiratory stability simultaneously:

```
all_conditions_check = (device_category == "IMV")
                     & (max_paralytics <= 0)
                     & (hemodynamic_stable == 1)
                     & (respiratory_stable == 1)
```

**Overnight quick-check flags:** `has_paralytics`, `has_hemo_stable`, `has_resp_stable`
**Failure UpSet columns:** `has_paralytics`, `hemo_unstable`, `resp_unstable`, `no_6h_window`

### 6.5 Summary Comparison Table

| Variant | `all_conditions_check` components | Extra null fills | Failure UpSet columns |
|---------|----------------------------------|------------------|-----------------------|
| Standard | IMV + no paralytics | None | 2 (`has_paralytics`, `no_6h_window`) |
| Hemodynamic | + NEE ≤ 0.2 | `nee` → 0.0 | 3 (+ `hemo_unstable`) |
| Respiratory | + FiO2 ≤ 0.5, PEEP ≤ 8, SpO2 ≥ 88 | `fio2_set` → 0.21, `peep_set` → 0.0, `spo2` → 100.0 | 3 (+ `resp_unstable`) |
| Both | + all of the above | All of the above | 4 (+ `hemo_unstable`, `resp_unstable`) |

---

## 7. Eligible Event Flagging

After the 6-hour check identifies eligible days, two flags are assigned to rows:

### 7.1 `eligible_event`

Set to `1.0` on exactly **one row per eligible day**: the first observation at or
after `event_time_at_6_hours`. This marks the moment eligibility was established and
SBT delivery assessment can begin.

### 7.2 Last-Observation Guard

If the `eligible_event` row happens to be the very last observation for that patient
(across all days), the flag is removed (set to null).

> **Why?** SBT delivery detection requires looking at a forward window (2 min or
> 30 min). If there is no subsequent data after the eligible event, we cannot
> determine whether a FLIP actually occurred and was sustained. Flagging such a row
> as eligible would produce a guaranteed false negative.

### 7.3 `on_vent_eligible`

Set to `1` on **all rows** belonging to days that have an `eligible_event`. This is
a convenience flag used to filter down to eligible-day rows for the delivery detection
loop.

---

## 8. SBT Delivery Detection

Delivery detection runs only on rows from eligible days (`on_vent_eligible == 1`).

### 8.1 FLIP Definition

A "FLIP" is a ventilator mode transition to a less supportive mode. The FLIP condition
is:

```
flip_check_flag = (device_category == "IMV")
                & (
                    (
                        (mode_category contains "pressure support" OR "cpap")
                        & (pressure_support_set <= 8)
                        & (peep_set <= 8)
                    )
                    OR (mode_name matches "^t-?piece$")
                  )
```

Key details:
- `mode_category` and `mode_name` are checked case-insensitively
- Null `pressure_support_set` and `peep_set` are filled with `inf` (treated as failing)
- T-piece detection uses a regex pattern matching `t-piece` or `tpiece`
- The FLIP must occur **after** the 6-hour eligibility threshold time (`_elig_threshold_time`)

### 8.2 Two Sustained Checks

At each candidate FLIP timepoint (where `flip_check_flag` is true and time > threshold),
the algorithm looks forward and checks whether the FLIP is sustained:

| Flag | Forward window | Condition |
|------|---------------|-----------|
| `SBT_delivery_2min` | 2 minutes | All observations within the window have `flip_check_flag == True` |
| `SBT_delivery_30min` | 30 minutes | All observations within the window have `flip_check_flag == True` |

### 8.3 One Detection Per Day

The algorithm finds at most one 2-min delivery and one 30-min delivery per day:
- Once a 2-min FLIP is found, it also checks if that same point satisfies the 30-min
  window
- After detecting both flags (or exhausting candidates), the loop moves to the next day
- A value of `1.0` means delivery was detected; null means it was not

---

## 9. Day-Level Aggregation

The row-level flags are collapsed into one row per ventilator-day.

### 9.1 Vent-Day Filter

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

### 9.2 Key Aggregated Columns

| Column | Logic |
|--------|-------|
| `is_eligible` | `True` if `eligible_event == 1` on any row that day |
| `SBT_delivery_2min` | `max()` of the flag across all rows in the day (1 if any row flagged) |
| `SBT_delivery_30min` | Same as above for 30-min delivery |
| `sbt_ground_truth` | `1.0` if **any** of the 4 EHR flowsheet columns is 1 **and** the day was eligible; null otherwise |

### 9.3 Ground Truth Construction

The ground truth is constructed from four EHR flowsheet columns:

```
sbt_ground_truth = 1.0 when (
    sbt_screen_pass_fail == 1
    OR sbt_screen_performed == 1
    OR sbt_delivery_pass_fail == 1
    OR sbt_delivery_performed == 1
) AND eligible_event == 1
```

If none of these conditions are met, `sbt_ground_truth` is null.

### 9.4 Delivery Failure Reasons

For eligible days where delivery was not detected:

| Column | Meaning |
|--------|---------|
| `SBT_delivery_2min_failure` | `"No mode FLIP to PS/CPAP/T-piece"` if no FLIP occurred at all; otherwise `"FLIP not sustained for 2min"` |
| `SBT_delivery_30min_failure` | `"No mode FLIP to PS/CPAP/T-piece"` if no FLIP occurred at all; otherwise `"FLIP not sustained for 30min"` |

---

## 10. Concordance Analysis

The algorithmic flags are compared against the EHR ground truth using standard
classification metrics.

### 10.1 Metrics Computed

For each of `SBT_delivery_2min` and `SBT_delivery_30min` vs `sbt_ground_truth`:

- **Confusion matrix** (TP, FP, FN, TN)
- **Accuracy**, **Precision**, **Recall**, **F1**, **Specificity**
- **Cohen's Kappa** with 95% bootstrap confidence interval (2000 resamples, seed=42)

### 10.2 Landis-Koch Interpretation Scale

| Kappa range | Interpretation |
|-------------|---------------|
| < 0.00 | Poor |
| 0.00 - 0.20 | Slight |
| 0.21 - 0.40 | Fair |
| 0.41 - 0.60 | Moderate |
| 0.61 - 0.80 | Substantial |
| 0.81 - 1.00 | Almost Perfect |

### 10.3 Outputs

- `delivery_concordance_summary.csv` — all metrics in one table
- `confusion_matrix_{flag_name}.png` — JAMA-style heatmap per flag (cividis colormap)

---

## 11. UpSet / Failure Analysis

### 11.1 Eligibility Failures UpSet

For each variant, an UpSet plot visualizes the intersection of failure reasons among
non-eligible days that had overnight data and IMV. The boolean columns vary by variant:

| Variant | UpSet columns |
|---------|--------------|
| Standard | `has_paralytics`, `no_6h_window` |
| Hemodynamic | `has_paralytics`, `hemo_unstable`, `no_6h_window` |
| Respiratory | `has_paralytics`, `resp_unstable`, `no_6h_window` |
| Both | `has_paralytics`, `hemo_unstable`, `resp_unstable`, `no_6h_window` |

### 11.2 Detection Failures UpSet

For eligible days where the 30-min SBT was not detected, three boolean columns capture
the failure mode:

| Column | Definition |
|--------|-----------|
| `no_flip` | `SBT_delivery_2min_failure == "No mode FLIP to PS/CPAP/T-piece"` |
| `flip_not_sustained_2min` | `SBT_delivery_2min_failure == "FLIP not sustained for 2min"` |
| `flip_sustained_2min_not_30min` | `SBT_delivery_2min` is not null but `SBT_delivery_30min` is null (FLIP sustained 2 min but not 30 min) |

Both UpSet plots are produced using `plot_upset()` from `helper.py`, which saves
a PNG visualization and a companion CSV with intersection counts.

---

## 12. Output Files

Each variant writes to its own output subdirectories under `output_phi/` (PHI) and
`output_to_share/` (aggregate, no PHI).

| File | Location | Description |
|------|----------|-------------|
| `sbt_day_level.parquet` | `output_phi/{variant}/` | Day-level dataset: one row per vent-day with eligibility, delivery flags, ground truth, and failure reasons |
| `consort_sbt_{SITE}.json` | `output_to_share/{variant}/` | CONSORT attrition funnel as structured JSON |
| `consort_sbt_{SITE}.png` | `output_to_share/{variant}/` | CONSORT diagram visualization |
| `delivery_concordance_summary.csv` | `output_to_share/{variant}/` | Concordance metrics (kappa, accuracy, precision, recall, F1, etc.) |
| `confusion_matrix_SBT_delivery_2min.png` | `output_to_share/{variant}/` | Confusion matrix heatmap for 2-min delivery |
| `confusion_matrix_SBT_delivery_30min.png` | `output_to_share/{variant}/` | Confusion matrix heatmap for 30-min delivery |
| `upset_eligibility_{SITE}.png` | `output_to_share/{variant}/` | UpSet plot of eligibility failure intersections |
| `upset_eligibility_{SITE}.csv` | `output_to_share/{variant}/` | Eligibility failure intersection counts |
| `upset_detection_{SITE}.png` | `output_to_share/{variant}/` | UpSet plot of detection failure intersections |
| `upset_detection_{SITE}.csv` | `output_to_share/{variant}/` | Detection failure intersection counts |
| `eligibility_failures_{SITE}.csv` | `output_phi/{variant}/` | Per-day eligibility failure flags (PHI) |
| `detection_failures_{SITE}.csv` | `output_phi/{variant}/` | Per-day detection failure flags (PHI) |
| `sbt_table1_{SITE}.csv` | `output_to_share/{variant}/` | Table 1 summary statistics |
| `sbt_table1_{SITE}.json` | `output_to_share/{variant}/` | Table 1 machine-readable summary |

**Variant subdirectories:**
- `sbt_standard/`
- `sbt_hemodynamic/`
- `sbt_respiratory/`
- `sbt_both_stabilities/`

---

## 13. Data Flow Diagram

```
                        ┌───────────────────────────┐
                        │   wide_dataset.parquet     │
                        │  (from 02_wide_dataset.py) │
                        └─────────────┬─────────────┘
                                      │
                        ┌─────────────v─────────────────┐
                        │  Stability Flag Computation    │
                        │  (variant-dependent)           │
                        │  hemodynamic_stable, resp...   │
                        └─────────────┬─────────────────┘
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
                    │   Quick-check: IMV? No paralytics?   │
                    │   (+stability flags per variant)     │
                    │   Hierarchical failure reasons        │
                    │   + CONSORT counts                    │
                    └──────────────┬──────────────────────┘
                                   │
                        passes screening
                                   │
                                   v
                ┌──────────────────────────────────────────┐
                │  6-Hour Continuous Eligibility Check       │
                │  Contiguous segments where:                │
                │    IMV + no paralytics                     │
                │    (+ stability conditions per variant)    │
                │  First segment >= 6h → eligible            │
                │  Records event_time_at_6_hours             │
                └──────────────┬───────────────────────────┘
                               │
                        eligible days
                               │
                               v
                ┌──────────────────────────────────────┐
                │   Eligible Event Flagging              │
                │   eligible_event = first obs >= 6h     │
                │   mark on that row (with last-obs      │
                │   guard), on_vent_eligible on all      │
                │   rows of eligible days                │
                └──────────────┬───────────────────────┘
                               │
                               v
                ┌──────────────────────────────────────────┐
                │   SBT Delivery Detection (FLIP)           │
                │   Mode transition to PS/CPAP/T-piece      │
                │   with PS ≤ 8 and PEEP ≤ 8                │
                │   Two checks:                              │
                │     SBT_delivery_2min  (sustained 2 min)   │
                │     SBT_delivery_30min (sustained 30 min)  │
                └──────────────┬───────────────────────────┘
                               │
                               v
                ┌──────────────────────────────────────────┐
                │   Day-Level Aggregation                   │
                │   Vent-day filter (post-intubation +      │
                │   cross-day IMV join)                      │
                │   is_eligible, delivery flags (max),       │
                │   sbt_ground_truth, failure reasons         │
                └──────┬──────────────┬───────────────────┘
                       │              │
              ┌────────┘              └──────────┐
              v                                  v
   ┌────────────────────┐           ┌───────────────────────┐
   │  Concordance        │           │  UpSet / Failure      │
   │  Cohen's Kappa +    │           │  Analysis             │
   │  bootstrap CI,      │           │  Eligibility +        │
   │  confusion matrices │           │  Detection failures   │
   └────────┬───────────┘           └───────────┬───────────┘
            │                                    │
            v                                    v
   delivery_concordance_summary.csv     upset_*.png / .csv
   confusion_matrix_*.png

              ┌─────────────────────────────────┐
              │   Final Outputs                  │
              │   sbt_day_level.parquet           │
              │   consort_sbt_{SITE}.json / .png  │
              │   sbt_table1_{SITE}.csv / .json   │
              └─────────────────────────────────┘
```

---

## 14. Key Design Decisions

### Why 6 hours, not 4 hours like SAT?

SBT eligibility uses a 6-hour continuous window (vs SAT's 4-hour window). The longer
window reflects the clinical rationale that ventilator weaning readiness requires a
more sustained period of stability on the ventilator without paralysis, ensuring the
patient has been stable on IMV long enough to attempt spontaneous breathing.

### Why no sedation requirement?

Unlike SAT (which requires `min_sedation_dose_2 > 0` — active sedation), SBT does
not check for sedation at all. This is because SBT assesses ventilator weaning
readiness, not sedation status. A patient can be eligible for an SBT regardless of
whether they are receiving sedation, as the trial tests respiratory capability rather
than neurological awakening.

### Why does the FLIP definition use low PS/PEEP thresholds?

The FLIP requires `pressure_support_set <= 8` and `peep_set <= 8`. These thresholds
define "minimal support" — the ventilator is providing just enough to offset the
resistance of the endotracheal tube while the patient does most of the breathing work.
Higher PS or PEEP would indicate the ventilator is still providing substantial
assistance, which would not constitute a true spontaneous breathing trial.

### Why two delivery thresholds (2 min and 30 min)?

- **2-minute threshold**: captures any mode FLIP that was sustained briefly. This is a
  sensitive detector that catches even short FLIPs that may represent SBT attempts that
  were quickly aborted.
- **30-minute threshold**: captures sustained FLIPs that represent a clinically
  meaningful SBT. Most institutional protocols define an SBT as lasting 30 minutes
  to 2 hours.

Having both allows analysis of the concordance gap between brief and sustained mode
changes, and lets researchers choose the threshold appropriate for their definition.

### Why exclude SAT/SBT flowsheet columns from forward-fill?

These are point-in-time charting events. Forward-filling would propagate a single
nurse documentation event across all subsequent timestamps, inflating the ground truth.

### Why is the cross-day IMV join needed?

Calendar-day grouping assigns a 22:30 observation to the previous calendar day's
`hosp_id_icu_day`. But that observation falls within the current day's 22:00-06:00
overnight window. The cross-day join re-checks IMV using clock-time boundaries rather
than grouping keys.

### Why only the first qualifying 6-hour segment per day?

One eligible segment is sufficient to establish SBT eligibility. Taking the earliest
one maximizes the remaining time window available for delivery detection.

### Why the last-observation guard on `eligible_event`?

Delivery detection requires a forward window (up to 30 minutes). If the eligible event
is the patient's very last observation, there is no data to evaluate delivery against —
the result would be an automatic false negative. Removing the flag prevents this
artifact.
