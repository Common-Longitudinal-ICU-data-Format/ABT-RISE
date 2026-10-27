# ETT-IMV proxy comparison

## Interpretation

This is an electronic-evidence coverage and robustness comparison, not a validation of clinical intubation/extubation accuracy.
The current code remains unchanged. Candidate timestamps are EHR-derived proxies, not verified ETT insertion/removal times.
Counts below 11 are withheld and other counts are rounded to the nearest 10 for this public diagnostic report. Exact aggregates and records remain private. This conservative display rule does not replace institutional disclosure review.

## Fixed baseline study panel

| Method | Confirmed starts | Confirmed ends | Unconfirmed ends | Zero elapsed intervals | Intervals with gaps | Median start shift, min [IQR] | Median end shift, min [IQR] |
|---|---:|---:|---:|---:|---:|---:|---:|
| B0 | ~9,740 | ~8,430 | ~1,300 | ~1,450 | — | 0 [0, 0] | 0 [0, 0] |
| A0 | ~9,740 | ~8,010 | ~1,720 | ~30 | — | 0 [0, 0] | 0 [0, 41.983] |
| C1 | ~9,440 | ~5,910 | ~3,530 | 0 | ~5,570 | 0 [0, 0] | 15 [0.017, 60] |
| C2 | ~9,700 | ~6,230 | ~3,460 | 0 | ~5,920 | 0 [0, 0] | 15 [0.017, 60] |
| C3 | ~9,700 | ~6,170 | ~3,530 | 0 | ~5,950 | 0 [0, 0] | 15 [0.017, 60] |
| C4 | ~9,700 | ~6,240 | ~3,460 | 0 | ~5,930 | 0 [0, 0] | 15 [0.017, 60] |

The denominator is the existing algorithm-defined study cohort; retention is not clinical sensitivity. B0/A0 use the baseline waterfall, whereas C1-C4 use original observations after the existing loader's timestamp conversion. Cross-family differences combine source and rule differences; they are not attributed solely to persistence thresholds.
Timing diagnostics use exact elapsed minutes. The baseline's stored whole-hour duration field is unchanged.

## Definitions

- **B0:** Unchanged current code.
- **A0:** Baseline with oxygen-flow input neutralized.
- **C1:** Two observations over 30 min for start and end.
- **C2:** First explicit IMV; two off observations over 30 min.
- **C3:** First explicit IMV; two off observations over 60 min.
- **C4:** C2 plus separately tagged probable IMV evidence.

C4 requires at least two probable observations spanning 30 minutes for inferred onset; it never overrides an explicitly non-IMV device. PEEP alone and PS/CPAP mode alone do not establish IMV.

## Additional discovery

The non-respiratory eligibility rules are inherited from the original pre-event ICU panel. Only cases outside the baseline cohort are counted here; additional evidence is not proof of additional correct cases.

| Method | New possible admissions | Probable-onset admissions |
|---|---:|---:|
| C1 | 0 | 0 |
| C2 | 0 | 0 |
| C3 | 0 | 0 |
| C4 | ~480 | ~480 |

## Data quality and timing limitations

- Maximum evidence gap for the main comparison: 120 minutes; this is an exploratory operational parameter, not a validated clinical standard.
- Median interval between known-state observations: 60 [49, 79] minutes; hourly observations can make 30/60-minute confirmation contrasts coincide.
- Source timezone interpretation remains unvalidated. Existing loader conversions and the baseline's removal of timezone metadata are preserved; these results are exploratory.
- A documentation gap is not an observed extubation. Private segment metadata retains uncertainty; the cohort export uses the first confirmed onset and earliest confirmed cessation after it in the selected window, without asserting continuous exposure through gaps.
- Ventilation via tracheostomy is not ETT ventilation. New tracheostomy censors ETT-specific eligibility without creating an extubation event; unknown route is not assumed absent.
- A mode change to pressure support/CPAP can remain IMV. T-piece/unsupported intervals can reflect breathing trials; no existing SBT classification was used to construct boundaries.
- Baseline method-label ties can vary while their event timestamps agree. The replay accepts a label difference only when original branch candidates prove equal timestamps.
- Routine categories may themselves be ETL-derived; ignoring dedicated event columns does not establish that upstream mappings were independent of human documentation.

## Scientific labels

Recommended primary terminology: **rule-based EHR-derived signatures compatible with SAT and SBT delivery**. For the cohort: **an EHR-derived endotracheal-IMV phenotype without recorded tracheostomy**. For timing: **EHR-derived proxies for IMV initiation and discontinuation**.
Do not claim exclusively device/pump-origin data, verified airway events, measured clinical recall, or superiority to nurse documentation without lineage verification and an independent validation design.

## Separate SAT/SBT handoff

Private cohort exports preserve the incumbent schema and add proxy confirmation, interval bounds, gap flags, and eligible-end metadata. They are not a validated drop-in analysis replacement: regenerate the hourly/supporting datasets and explicitly honor the selected ETT eligibility bounds in a separate SAT/SBT run.
New discovery cases require all supporting datasets to be rebuilt. Existing SAT/SBT and R code was not changed or tuned to favor a proxy definition.

## Operational recommendation

Use the comparison to prioritize a high-coverage direct-evidence cohort and keep probable-only discovery separately identified. Assess coverage together with brief intervals, gaps, censoring, and timing shifts; do not choose by coverage alone or by favorable SAT/SBT associations. No uniquely most accurate definition can be selected without reference evidence.
