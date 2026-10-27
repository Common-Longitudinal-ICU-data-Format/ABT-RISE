# Scientific labels for ABT-RISE and its ventilation cohort

## Recommended terminology

Use **EHR-derived electronic signatures of spontaneous awakening and breathing trials** for the paper's main construct.
In Methods and results, use **SAT/SBT-delivery-compatible electronic signatures** until independent validation supports a stronger occurrence claim.

A suitable working title is:

> ABT-RISE: Rule-based electronic signatures of spontaneous awakening and breathing trials from structured ICU EHR data

For cohort discovery, use **an EHR-derived endotracheal-IMV phenotype without recorded tracheostomy**.
For timing, use **EHR-derived proxies for IMV initiation and discontinuation**, or **intubation/extubation timing proxies** with an explicit definition.
Do not describe a device-state transition as independently verified ETT insertion/removal.

“Direct evidence” means an explicit structured IMV device classification, not direct clinical observation or a validated physical airway event.
“Probable evidence” means a compatible ventilator pattern when the device classification is missing; it remains a distinct sensitivity/discovery class.

## Why these labels fit the clinical construct

Clinical guidance permits spontaneous breathing trials with or without pressure support. [1]
A randomized trial compared a 30-minute pressure-support trial with a two-hour T-piece trial before assessing successful extubation. [2]
A controlled-mode-to-pressure-support change can therefore be consistent with SBT initiation while the patient remains invasively ventilated; it is not itself extubation.

The current SBT implementation detects eligible-day transitions to low-assist PS/CPAP or T-piece conditions, with short stability windows.
The current SAT implementation detects eligible sedative/opioid or non-opioid infusion-interruption patterns over an observed forward window.
These are structured patterns compatible with trial delivery, not proof of clinician intent, completed safety screening, a completed trial, or successful liberation.
The existing algorithms are in `code/03_sat.py` and `code/04_sbt_both.py`; this cohort comparison did not change or retune them.

EHR initiation-timing validation has used institution-specific ventilator data and known invasive-airway context; PEEP alone was explicitly recognized as present during NIV too. [3]
That supports source-specific validation rather than treating a single setting or a 30-minute observation span as a universal clinical standard.
CLIF mCIDE standardizes device/mode categories and preserves source names, but neither the schema nor a category name certifies direct machine provenance. [4]

## What the local comparison supports

The unchanged current cohort code was run as a byte-identical copy in private staging.
Its replay matched onset/cessation timestamps, durations, and no-end reasons for the full baseline panel.
Same-time competing end branches can produce nondeterministic method labels; those label differences were accepted only when the original branch candidates proved identical times.

The public diagnostic counts below are rounded; exact aggregates remain private.
They describe electronic evidence under the existing study eligibility rules, not measured clinical sensitivity or specificity.

| Definition | Baseline-panel starts | Confirmed end proxies | Operational interpretation |
|---|---:|---:|---|
| Current code, unchanged | ~9,740 | ~8,430 | Incumbent comparator, including the oxygen-flow branch |
| Baseline oxygen-flow ablation | ~9,740 | ~8,010 | Isolates that branch on the incumbent input source |
| Two observations spanning 30 min for both states | ~9,440 | ~5,910 | More restrictive start confirmation |
| First explicit IMV; two off observations spanning 30 min | ~9,700 | ~6,230 | High-coverage direct-evidence entry candidate |
| First explicit IMV; two off observations spanning 60 min | ~9,700 | ~6,170 | End-persistence sensitivity, not demonstrated greater accuracy |
| Evidence-assisted onset; 30-min off confirmation | ~9,700 | ~6,240 | Adds tagged probable evidence; finds ~480 additional candidate admissions outside the baseline panel |

Immediate explicit-IMV entry retains approximately **99.6%** of the baseline's algorithm-defined population.
This is **baseline retention**, not evidence that 99.6% of true clinical cases were found.
The oxygen-flow ablation greatly reduces zero-elapsed intervals, but eliminating that pattern alone does not validate clinical event times.

### Documentation-gap sensitivity is material

| Maximum gap, min | Direct-entry starts | Confirmed ends with 30-min off persistence |
|---|---:|---:|
| 60 | ~9,700 | ~5,020 |
| 120 | ~9,700 | ~6,230 |
| 240 | ~9,700 | ~7,090 |

Entry is robust across these tested gap assumptions; end availability is not.
The median interval between known-state observations is approximately 60 minutes, so 30- and 60-minute confirmation can occur at the same next observation.
A larger gap allowance yields more end proxies but does not prove the intervening state was continuous.

At the 120-minute gap setting, only about **64%** of retained direct-entry cases have a confirmed end proxy within the analysis window.
Do not exclude the remainder from cohort membership just because end confirmation is missing.
Retain right-censored/uncertain endings and reported interval bounds; missing confirmation does not mean extubation never occurred.
The code's onset-left-censor flag means no prior known non-IMV observation was found in the window; it does not automatically prove the patient was intubated before hospital arrival.

## Provisional operational recommendation

1. **Primary discovery population:** prioritize directly documented IMV entry with recorded tracheostomy absence rather than requiring 30 minutes of start confirmation for every case.
2. **Sensitivity/discovery extension:** retain ventilator-pattern-only probable cases separately. These additional candidates warrant validation; do not silently label them confirmed ETT cases.
3. **End timing:** use the sustained known non-IMV transition as a declared proxy, retain uncertain bounds, and report the gap/persistence sensitivities. This run does not identify a uniquely most accurate end definition.
4. **Do not rank by coverage alone:** examine brief intervals, documentation gaps, censoring, observed resumptions and source contradictions together. Do not choose a definition because it produces a desired SAT/SBT association.

This recommendation favors high-coverage cohort entry while making weaker route/timing evidence visible.
It is not authorization to replace the production definition or a declaration that the new method is clinically superior.

## Separate SAT/SBT validation

Keep cohort construction and SAT/SBT detection separate to avoid circular validation.
The proxy detector must not use SAT/SBT classifications to decide ventilation boundaries, and the signature algorithms must not be tuned to make one cohort definition look better.

Recommended construct labels are:

- **Electronic eligibility phenotype:** a rule-defined period in which a trial could be evaluated; not evidence that a trial happened.
- **SAT-compatible infusion-interruption signature:** an observed medication pattern satisfying the declared SAT rule.
- **SBT-initiation-compatible respiratory signature:** a transition satisfying the declared low-assist/spontaneous-breathing rule; not proof of a completed or successful SBT.
- **Nurse-charted trial documentation:** the human criterion comparator, not automatically a perfect clinical gold standard.
- **Independent clinical/device reference:** separately adjudicated evidence used when claiming event accuracy or superiority to documentation.

An electronic-positive/chart-negative case could represent a true undocumented trial, a nonspecific interruption or mode change, or bad source data.
Do not call it a clinical false positive solely because the chart lacks a trial label.
Likewise, structured data may improve capture of actual intervention changes, but it cannot by itself establish clinician intent.

The existing code uses “ground truth” in some comparison comments; the scientific paper should instead describe **agreement with nurse-charted documentation** unless the reference has been independently adjudicated.
Before claiming machine-derived signals or superiority, verify field lineage, timestamp provenance, coverage of confirmation windows, and the distinction between missing values and true zero/absence.

## Downstream use and limitations

Private variant exports preserve the cohort schema and include confirmation status, interval bounds, evidence tags, gap flags and `proxy_eligible_end`.
Use those with regenerated hourly/supporting data in a separate SAT/SBT run; do not simply swap a cohort file into baseline intermediates.
New discovery cases need their supporting clinical datasets rebuilt.
Existing duration/outcome calculations also need a censoring-aware review before missing proxy ends are interpreted as no clinical extubation or assigned a ventilation-free-day outcome.

No clinical reference labels were available in this run.
Source timezone interpretation remains unvalidated; baseline loader conversions and its timezone-stripping behavior were preserved.
Structured `_category` fields may already incorporate upstream human documentation or derived flags, so computational independence from event columns does not establish full ETL independence.
There is no demonstrated SAT/SBT equivalence/superiority result from this cohort comparison.

## Reproducible local workflow

All clinical records, logs, exact aggregates, and variant cohort exports remain under the Git-ignored `output_phi/imv_comparison/` tree.
Do not upload those files to third-party review tools or public repositories.
The generated public diagnostic report is `docs/reports/ett-imv-proxy-comparison.md`; its rounding/suppression is a conservative display policy, not a substitute for institutional disclosure review.

```bash
uv sync --locked
uv run python -m unittest discover -s tests -v
uv run python code/imv_baseline.py
uv run python code/compare_imv_definitions.py --gap-sensitivity
```

The baseline runner supports an explicitly supplied existing waterfall cache and records that source; it never rewrites the incumbent notebook.
The comparison rechecks source/event parity and protects input files against output-path collisions.
The original `code/01_cohort.py`, SAT/SBT modules, and R analysis scripts were not changed by this work.

## References

[1] Roberts KJ et al. “AARC Clinical Practice Guideline: Spontaneous Breathing Trials for Liberation From Adult Mechanical Ventilation.” *Respiratory Care* 69, 891–901 (2024). doi:10.4187/respcare.11735. Abstract reviewed.
https://paperclip.gxl.ai/citations/papers/154830580#L1

[2] Subirà C et al. “Effect of Pressure Support vs T-Piece Ventilation Strategies During Spontaneous Breathing Trials on Successful Extubation Among Patients Receiving Mechanical Ventilation: A Randomized Clinical Trial.” *JAMA* 321, 2175–2182 (2019). doi:10.1001/jama.2019.7234. Abstract reviewed.
https://paperclip.gxl.ai/citations/papers/31471125#L1

[3] Smischney NJ et al. “Derivation and validation of a search algorithm to retrospectively identify mechanical ventilation initiation in the intensive care unit.” *BMC Medical Informatics and Decision Making* 14, 55 (2014). doi:10.1186/1472-6947-14-55.
https://paperclip.gxl.ai/citations/papers/PMC4082421#L27-L44

[4] CLIF Consortium. Respiratory-support dictionary and mCIDE device/mode specifications, CLIF 2.1.
https://clif-icu.com/data-dictionary/data-dictionary-2.1.0
