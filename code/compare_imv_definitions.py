"""Compare the unchanged baseline with independent respiratory EHR proxies.

Clinical records, exact aggregates, and cohort exports remain private. The public
Markdown report is rounded/redacted and is not a clinical-accuracy evaluation.
"""
from __future__ import annotations

import argparse
from dataclasses import replace, asdict
import json
from pathlib import Path
import sys
import tempfile
import traceback
from typing import Any

import polars as pl

from imv_baseline import EVENT_FIELDS, assert_baseline_parity, file_sha256, private_dir, private_json, private_parquet, replay_baseline
from imv_proxy import DEFINITIONS, ProxyDefinition, detect_normalized_episodes, normalize_observations

KEYS = ["hospitalization_id", "icu_start", "window_start", "window_end"]
LABELS = {
    "B0": "Unchanged current code", "A0": "Baseline with oxygen-flow input neutralized",
    "C1": "Two observations over 30 min for start and end", "C2": "First explicit IMV; two off observations over 30 min",
    "C3": "First explicit IMV; two off observations over 60 min", "C4": "C2 plus separately tagged probable IMV evidence",
}


def aggregate_proxy_boundaries(episodes: pl.DataFrame, tracheostomy_points: pl.DataFrame | None = None) -> pl.DataFrame:
    """First onset/earliest cessation per eligible window, with gap uncertainty.

A data gap is not an affirmative extubation. A later observed cessation can still
bound the first onset; that interval is not asserted to be continuous exposure.
"""
    confirmed = episodes.filter(pl.col("start_confirmed"))
    eligible_keys = confirmed.filter(pl.col("has_icu_imv_evidence")).select(KEYS).unique()
    confirmed = confirmed.join(eligible_keys, on=KEYS, how="semi").sort(*KEYS, "start_time")
    starts = confirmed.group_by(KEYS, maintain_order=True).agg(
        pl.col("start_time").first(), pl.col("start_confirmation_time").first(),
        pl.col("start_evidence").first(), pl.col("start_lower_bound").first(), pl.col("left_censored").first(),
        pl.len().alias("segment_count"),
        pl.when(pl.col("censor_reason") == "new_tracheostomy").then(pl.col("censor_time")).otherwise(None).min().alias("trach_censor_time"),
    )
    ends = (confirmed.filter(pl.col("end_time").is_not_null()).sort(*KEYS, "end_time")
            .unique(subset=KEYS, keep="first")
            .select(*KEYS, "end_time", "end_confirmation_time", "end_lower_bound", "end_upper_bound"))
    result = starts.join(ends, on=KEYS, how="left")
    if tracheostomy_points is not None:
        points = (tracheostomy_points.filter(pl.col("_trach") == 1)
                  .select("hospitalization_id", "recorded_dttm")
                  .join(eligible_keys, on="hospitalization_id", how="inner")
                  .filter(pl.col("recorded_dttm").is_between(pl.col("window_start"), pl.col("window_end")))
                  .group_by(KEYS).agg(pl.col("recorded_dttm").min().alias("source_trach_time")))
        result = result.join(points, on=KEYS, how="left").with_columns(
            pl.min_horizontal("trach_censor_time", "source_trach_time").alias("trach_censor_time")
        ).drop("source_trach_time")
    interval_rows = (episodes.join(result.select(*KEYS, pl.col("start_time").alias("first_start"), pl.col("end_time").alias("first_end")), on=KEYS, how="inner")
                     .filter((pl.col("start_time") >= pl.col("first_start"))
                             & (pl.col("start_time") <= pl.coalesce("first_end", "window_end"))))
    flags = interval_rows.group_by(KEYS).agg(
        (pl.col("censor_reason") == "documentation_gap").any().alias("interval_has_gap"),
        pl.col("had_unknown_observations").any().alias("interval_has_unknown"),
        pl.col("had_conflicting_observations").any().alias("interval_has_conflict"),
        pl.col("used_probable_evidence").any().alias("used_probable_evidence"),
        pl.col("end_upper_bound").drop_nulls().min().alias("observed_off_time"),
    )
    return result.join(flags, on=KEYS, how="left").sort(*KEYS)


def _export_cohort(cohort: pl.DataFrame, boundaries: pl.DataFrame, definition: str) -> pl.DataFrame:
    """Keep the incumbent schema and add explicit proxy status/eligible boundaries."""
    keys = cohort.select("hospitalization_id", pl.col("first_icu_start").alias("icu_start"),
                         pl.col("admission_dttm").alias("window_start"), pl.col("first_icu_end").alias("window_end"))
    keys = keys.with_columns(*[pl.col(name).cast(boundaries.schema[name]) for name in KEYS if name != "hospitalization_id"])
    selected = keys.join(boundaries, on=KEYS, how="left")
    replacements = selected.select(
        "hospitalization_id", pl.col("start_time").alias("intubation_time"), pl.col("end_time").alias("extubation_time"),
        pl.col("start_confirmation_time").alias("proxy_start_confirmation_time"),
        pl.col("end_confirmation_time").alias("proxy_end_confirmation_time"),
        pl.col("start_lower_bound").alias("proxy_start_lower_bound"),
        pl.col("end_lower_bound").alias("proxy_end_lower_bound"), pl.col("end_upper_bound").alias("proxy_end_upper_bound"),
        pl.col("start_evidence").alias("proxy_start_evidence"), pl.col("left_censored").alias("proxy_start_left_censored"),
        pl.col("trach_censor_time").alias("proxy_tracheostomy_censor_time"),
        pl.col("interval_has_gap").alias("proxy_interval_has_gap"), pl.col("interval_has_unknown").alias("proxy_interval_has_unknown"),
        pl.col("interval_has_conflict").alias("proxy_interval_has_conflict"), pl.col("segment_count").alias("proxy_segment_count"),
        pl.col("observed_off_time").alias("proxy_observed_off_time"),
    )
    base = cohort.drop([c for c in ("intubation_time", "extubation_time", "imv_duration_hours", "extubation_method", "no_extubation_reason", "device_before_intubation", "device_after_extubation") if c in cohort.columns])
    exported = base.join(replacements, on="hospitalization_id", how="left").with_columns(
        pl.lit(definition).alias("proxy_definition"), pl.col("intubation_time").is_not_null().alias("proxy_has_confirmed_start"),
        pl.col("extubation_time").is_not_null().alias("proxy_has_confirmed_end"),
        pl.min_horizontal("first_icu_end", "proxy_tracheostomy_censor_time", "extubation_time").alias("proxy_eligible_end"),
        ((pl.col("extubation_time") - pl.col("intubation_time")).dt.total_microseconds() / 3_600_000_000).alias("imv_duration_hours"),
        pl.when(pl.col("extubation_time").is_not_null()).then(pl.lit(definition)).otherwise(None).alias("extubation_method"),
        pl.when(pl.col("extubation_time").is_not_null()).then(pl.lit(None, dtype=pl.String))
        .when(pl.col("intubation_time").is_null()).then(pl.lit("No confirmed ETT-IMV start under this proxy"))
        .when(pl.col("proxy_tracheostomy_censor_time").is_not_null()).then(pl.lit("ETT-specific follow-up censored at new tracheostomy"))
        .otherwise(pl.lit("Unconfirmed IMV cessation within observation window")).alias("no_extubation_reason"),
        pl.lit(True).alias("proxy_duration_is_elapsed_not_continuous"),
    )
    return exported


def _distribution(series: pl.Series) -> dict[str, Any]:
    values = series.drop_nulls()
    if not values.len():
        return {"n": 0, "median": None, "q25": None, "q75": None}
    return {"n": values.len(), "median": round(float(values.median()), 3),
            "q25": round(float(values.quantile(0.25)), 3), "q75": round(float(values.quantile(0.75)), 3)}


def _metrics(frame: pl.DataFrame, baseline: pl.DataFrame) -> dict[str, Any]:
    elapsed = (frame["extubation_time"] - frame["intubation_time"]).dt.total_microseconds() / 60_000_000
    starts = frame["intubation_time"].is_not_null().sum()
    ends = frame["extubation_time"].is_not_null().sum()
    paired = frame.select("hospitalization_id", "intubation_time", "extubation_time").join(
        baseline.select("hospitalization_id", "intubation_time", "extubation_time"), on="hospitalization_id", suffix="_baseline")
    result = {
        "panel_admissions": frame.height, "confirmed_starts": int(starts), "confirmed_ends": int(ends),
        "unconfirmed_ends": int(starts - ends), "zero_elapsed_intervals": int((elapsed == 0).fill_null(False).sum()),
        "brief_under_30_minutes": int(((elapsed >= 0) & (elapsed < 30)).fill_null(False).sum()),
        "elapsed_minutes": _distribution(elapsed),
    }
    for name in ("intubation_time", "extubation_time"):
        shift = (paired[name] - paired[name + "_baseline"]).dt.total_microseconds() / 60_000_000
        result[name + "_shift_minutes"] = _distribution(shift)
    for column, label in (("proxy_interval_has_gap", "intervals_with_gaps"), ("proxy_start_left_censored", "left_censored_starts"),
                          ("proxy_interval_has_unknown", "intervals_with_unknown"), ("proxy_interval_has_conflict", "intervals_with_conflicts")):
        result[label] = int(frame[column].fill_null(False).sum()) if column in frame.columns else None
    result["probable_starts"] = int((frame["proxy_start_evidence"] == "probable").fill_null(False).sum()) if "proxy_start_evidence" in frame.columns else None
    return result


def _prepare_baseline(metadata: dict[str, Any]) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, Path]:
    if metadata.get("parity", {}).get("mismatches") != 0 or not metadata.get("original_files_unchanged"):
        raise RuntimeError("Baseline validation is incomplete or failed")
    source = Path(metadata["staged_project"]) / "code" / "01_cohort.py"
    if file_sha256(source) != metadata["notebook_sha256"]:
        raise RuntimeError("Baseline source snapshot changed")
    cohort = pl.read_parquet(metadata["cohort"])
    wf = pl.read_parquet(metadata["waterfall"], columns=["hospitalization_id", "recorded_dttm", "device_category", "lpm_set"])
    dtype = wf.schema["recorded_dttm"]
    if isinstance(dtype, pl.Datetime) and dtype.time_zone is not None:
        wf = wf.with_columns(pl.col("recorded_dttm").dt.replace_time_zone(None))
    wf = (wf.join(cohort.select("hospitalization_id", "admission_dttm", "first_icu_end"), on="hospitalization_id", how="inner")
          .filter(pl.col("recorded_dttm").is_between(pl.col("admission_dttm"), pl.col("first_icu_end")))
          .sort("hospitalization_id", "recorded_dttm"))
    replay, ties = replay_baseline(wf, cohort.select("hospitalization_id", "discharge_category"), source,
                                   expected_sha256=metadata["notebook_sha256"], return_branch_candidates=True)
    assert_baseline_parity(cohort, replay, method_ties=ties)
    ablation_input = wf.with_columns(pl.when(pl.col("device_category") == "imv").then(None).otherwise(pl.col("lpm_set")).alias("lpm_set"))
    ablation = replay_baseline(ablation_input, cohort.select("hospitalization_id", "discharge_category"), source)
    ablation_cohort = (cohort.drop([c for c in EVENT_FIELDS if c != "hospitalization_id"])
                      .join(ablation.select(EVENT_FIELDS), on="hospitalization_id", how="left"))
    windows = cohort.select("hospitalization_id", pl.col("admission_dttm").alias("window_start"),
                            pl.col("first_icu_start").alias("icu_start"), pl.col("first_icu_end").alias("window_end"))
    return cohort, ablation_cohort, windows, source


def _quality(observations: pl.DataFrame, normalized: pl.DataFrame) -> dict[str, Any]:
    known = normalized.filter(pl.col("_signal") >= 0).with_columns(
        (pl.col("recorded_dttm").diff().over("hospitalization_id").dt.total_microseconds() / 60_000_000).alias("gap_minutes")
    )
    return {
        "original_rows": observations.height, "distinct_timestamp_rows": normalized.height,
        "input_missing_timestamp_rows": observations["recorded_dttm"].null_count(),
        "input_missing_key_rows": observations["hospitalization_id"].null_count(),
        "conflicting_device_timestamps": int(normalized["_conflict"].sum()),
        "tracheostomy_zero_rows": int((normalized["_trach"] == 0).fill_null(False).sum()),
        "tracheostomy_positive_rows": int((normalized["_trach"] == 1).fill_null(False).sum()),
        "tracheostomy_unknown_rows": normalized["_trach"].null_count(),
        "known_state_gap_minutes": _distribution(known.filter(pl.col("gap_minutes") > 0)["gap_minutes"]),
    }


def _public_count(number: int | None) -> str:
    if number is None:
        return "—"
    if 0 < number < 11:
        return "withheld"
    if number == 0:
        return "0"
    return f"~{int(round(number / 10) * 10):,}"


def _public_median(distribution: dict[str, Any]) -> str:
    if distribution["n"] < 11 or distribution["median"] is None:
        return "—"
    return f"{distribution['median']:g} [{distribution['q25']:g}, {distribution['q75']:g}]"


def render_report(summary: dict[str, Any]) -> str:
    lines = ["# ETT-IMV proxy comparison", "", "## Interpretation", "",
        "This is an electronic-evidence coverage and robustness comparison, not a validation of clinical intubation/extubation accuracy.",
        "The current code remains unchanged. Candidate timestamps are EHR-derived proxies, not verified ETT insertion/removal times.",
        "Counts below 11 are withheld and other counts are rounded to the nearest 10 for this public diagnostic report. Exact aggregates and records remain private. This conservative display rule does not replace institutional disclosure review.", "",
        "## Fixed baseline study panel", "",
        "| Method | Confirmed starts | Confirmed ends | Unconfirmed ends | Zero elapsed intervals | Intervals with gaps | Median start shift, min [IQR] | Median end shift, min [IQR] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for name, metric in summary["shared"].items():
        lines.append(f"| {name} | {_public_count(metric['confirmed_starts'])} | {_public_count(metric['confirmed_ends'])} | {_public_count(metric['unconfirmed_ends'])} | {_public_count(metric['zero_elapsed_intervals'])} | {_public_count(metric['intervals_with_gaps'])} | {_public_median(metric['intubation_time_shift_minutes'])} | {_public_median(metric['extubation_time_shift_minutes'])} |")
    lines += ["", "The denominator is the existing algorithm-defined study cohort; retention is not clinical sensitivity. B0/A0 use the baseline waterfall, whereas C1-C4 use original observations after the existing loader's timestamp conversion. Cross-family differences combine source and rule differences; they are not attributed solely to persistence thresholds.",
              "Timing diagnostics use exact elapsed minutes. The baseline's stored whole-hour duration field is unchanged.", "", "## Definitions", ""]
    lines += [f"- **{name}:** {label}." for name, label in LABELS.items() if name in summary["shared"]]
    lines += ["", "C4 requires at least two probable observations spanning 30 minutes for inferred onset; it never overrides an explicitly non-IMV device. PEEP alone and PS/CPAP mode alone do not establish IMV.", "",
              "## Additional discovery", "", "The non-respiratory eligibility rules are inherited from the original pre-event ICU panel. Only cases outside the baseline cohort are counted here; additional evidence is not proof of additional correct cases.", "",
              "| Method | New possible admissions | Probable-onset admissions |", "|---|---:|---:|"]
    for name, metric in summary["discovery"].items():
        lines.append(f"| {name} | {_public_count(metric['new_cases'])} | {_public_count(metric['new_probable_cases'])} |")
    lines += ["", "## Data quality and timing limitations", "",
              f"- Maximum evidence gap for the main comparison: {summary['max_gap_minutes']:g} minutes; this is an exploratory operational parameter, not a validated clinical standard.",
              f"- Median interval between known-state observations: {_public_median(summary['quality']['known_state_gap_minutes'])} minutes; hourly observations can make 30/60-minute confirmation contrasts coincide.",
              "- Source timezone interpretation remains unvalidated. Existing loader conversions and the baseline's removal of timezone metadata are preserved; these results are exploratory.",
              "- A documentation gap is not an observed extubation. Private segment metadata retains uncertainty; the cohort export uses the first confirmed onset and earliest confirmed cessation after it in the selected window, without asserting continuous exposure through gaps.",
              "- Ventilation via tracheostomy is not ETT ventilation. New tracheostomy censors ETT-specific eligibility without creating an extubation event; unknown route is not assumed absent.",
              "- A mode change to pressure support/CPAP can remain IMV. T-piece/unsupported intervals can reflect breathing trials; no existing SBT classification was used to construct boundaries.",
              "- Baseline method-label ties can vary while their event timestamps agree. The replay accepts a label difference only when original branch candidates prove equal timestamps.",
              "- Routine categories may themselves be ETL-derived; ignoring dedicated event columns does not establish that upstream mappings were independent of human documentation.", "",
              "## Scientific labels", "",
              "Recommended primary terminology: **rule-based EHR-derived signatures compatible with SAT and SBT delivery**. For the cohort: **an EHR-derived endotracheal-IMV phenotype without recorded tracheostomy**. For timing: **EHR-derived proxies for IMV initiation and discontinuation**.",
              "Do not claim exclusively device/pump-origin data, verified airway events, measured clinical recall, or superiority to nurse documentation without lineage verification and an independent validation design.", "",
              "## Separate SAT/SBT handoff", "",
              "Private cohort exports preserve the incumbent schema and add proxy confirmation, interval bounds, gap flags, and eligible-end metadata. They are not a validated drop-in analysis replacement: regenerate the hourly/supporting datasets and explicitly honor the selected ETT eligibility bounds in a separate SAT/SBT run.",
              "New discovery cases require all supporting datasets to be rebuilt. Existing SAT/SBT and R code was not changed or tuned to favor a proxy definition.", "",
              "## Operational recommendation", "",
              "Use the comparison to prioritize a high-coverage direct-evidence cohort and keep probable-only discovery separately identified. Assess coverage together with brief intervals, gaps, censoring, and timing shifts; do not choose by coverage alone or by favorable SAT/SBT associations. No uniquely most accurate definition can be selected without reference evidence.", ""]
    return "\n".join(lines)


def run_comparison(metadata_path: Path, output_root: Path, report_path: Path, *, max_gap_minutes: float = 120,
                   discovery: bool = True, gap_sensitivity: bool = False) -> dict[str, Any]:
    metadata_path, output_root, report_path = map(lambda p: Path(p).resolve(), (metadata_path, output_root, report_path))
    replace(DEFINITIONS["C1"], max_gap_minutes=max_gap_minutes)  # Validate before any writes.
    metadata = json.loads(metadata_path.read_text())
    staged = Path(metadata["staged_project"]).resolve()
    protected = {metadata_path, *(Path(metadata[key]).resolve() for key in ("cohort", "waterfall", "original_observations", "eligible_windows")),
                 staged / "code" / "01_cohort.py"}
    protected.update(path for path in (staged / "clif_config.json", staged / "code" / "helper.py") if path.is_file())
    planned = {report_path, output_root / "summary.json", output_root / "provenance.json"}
    planned.update(output_root / "cohorts" / f"{name}.parquet" for name in ("B0", "A0", "C1", "C2", "C3", "C4"))
    planned.update(output_root / directory / f"{name}.parquet" for directory in ("episodes", "discovery") for name in ("C1", "C2", "C3", "C4"))
    planned = {path.resolve() for path in planned}
    if protected & planned:
        raise ValueError("Comparison output would overwrite a protected input")
    before = {str(path): file_sha256(path) for path in protected}
    cohort, ablation, shared_windows, source = _prepare_baseline(metadata)
    private_dir(output_root)
    observations = pl.read_parquet(metadata["original_observations"])
    normalized = normalize_observations(observations)
    panel = pl.read_parquet(metadata["eligible_windows"])
    eligible = panel.select("hospitalization_id", pl.col("admission_dttm").alias("window_start"),
                            pl.col("icu_in_dttm").alias("icu_start"), pl.col("icu_out_dttm").alias("window_end"))
    summary: dict[str, Any] = {"shared": {}, "discovery": {}, "quality": _quality(observations, normalized),
                               "max_gap_minutes": max_gap_minutes, "baseline_parity": metadata["parity"], "gap_sensitivity": {}}
    for name, frame in (("B0", cohort), ("A0", ablation)):
        summary["shared"][name] = _metrics(frame, cohort)
        private_parquet(output_root / "cohorts" / f"{name}.parquet", frame)
    other_windows = eligible.filter(~pl.col("hospitalization_id").is_in(cohort["hospitalization_id"].unique().implode()))
    for name in ("C1", "C2", "C3", "C4"):
        definition = replace(DEFINITIONS[name], max_gap_minutes=max_gap_minutes)
        events = detect_normalized_episodes(normalized, shared_windows, definition)
        boundaries = aggregate_proxy_boundaries(events, normalized)
        exported = _export_cohort(cohort, boundaries, name)
        private_parquet(output_root / "episodes" / f"{name}.parquet", events)
        private_parquet(output_root / "cohorts" / f"{name}.parquet", exported)
        summary["shared"][name] = _metrics(exported, cohort)
        summary["discovery"][name] = {"new_cases": 0, "new_probable_cases": 0}
        if discovery and other_windows.height:
            signals = normalized.filter(pl.col("_signal").is_in([1, 2] if definition.evidence_assisted else [1]))
            candidate_windows = (signals.select("hospitalization_id", "recorded_dttm").join(other_windows, on="hospitalization_id", how="inner")
                                 .filter(pl.col("recorded_dttm").is_between(pl.col("icu_start"), pl.col("window_end")))
                                 .select(KEYS).unique())
            if candidate_windows.height:
                expanded_events = detect_normalized_episodes(normalized, candidate_windows, definition)
                expanded = (aggregate_proxy_boundaries(expanded_events, normalized).sort("hospitalization_id", "icu_start", "start_time")
                            .unique(subset=["hospitalization_id"], keep="first"))
                private_parquet(output_root / "discovery" / f"{name}.parquet", expanded)
                summary["discovery"][name] = {"new_cases": expanded.height, "new_probable_cases": int((expanded["start_evidence"] == "probable").sum())}
        print(f"{name}: confirmed starts={summary['shared'][name]['confirmed_starts']}; confirmed ends={summary['shared'][name]['confirmed_ends']}; new candidate admissions={summary['discovery'][name]['new_cases']}")
    if gap_sensitivity:
        for gap in (60.0, 120.0, 240.0):
            if gap == max_gap_minutes:
                summary["gap_sensitivity"][str(int(gap))] = {name: summary["shared"][name] for name in ("C1", "C2", "C3", "C4")}
                continue
            summary["gap_sensitivity"][str(int(gap))] = {}
            for name in ("C1", "C2", "C3", "C4"):
                events = detect_normalized_episodes(normalized, shared_windows, replace(DEFINITIONS[name], max_gap_minutes=gap))
                export = _export_cohort(cohort, aggregate_proxy_boundaries(events, normalized), name)
                summary["gap_sensitivity"][str(int(gap))][name] = _metrics(export, cohort)
    if any(file_sha256(path) != digest for path, digest in ((Path(path), digest) for path, digest in before.items())):
        raise RuntimeError("A protected comparison input changed during execution")
    private_json(output_root / "summary.json", summary)
    private_json(output_root / "provenance.json", {"baseline_metadata": metadata, "input_sha256": before,
                                                  "definitions": {name: asdict(replace(DEFINITIONS[name], max_gap_minutes=max_gap_minutes)) for name in ("C1", "C2", "C3", "C4")}})
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(render_report(summary))
    return summary


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=root / "output_phi/imv_comparison/latest_baseline.json")
    parser.add_argument("--output-root", type=Path, default=root / "output_phi/imv_comparison/definition_comparison")
    parser.add_argument("--report", type=Path, default=root / "docs/reports/ett-imv-proxy-comparison.md")
    parser.add_argument("--max-gap-minutes", type=float, default=120)
    parser.add_argument("--no-discovery", action="store_true")
    parser.add_argument("--gap-sensitivity", action="store_true")
    args = parser.parse_args()
    try:
        run_comparison(args.metadata, args.output_root, args.report, max_gap_minutes=args.max_gap_minutes,
                       discovery=not args.no_discovery, gap_sensitivity=args.gap_sensitivity)
    except Exception:
        try:
            directory = private_dir(args.output_root)
            with tempfile.NamedTemporaryFile(mode="w", prefix="comparison-error-", suffix=".log", dir=directory, delete=False) as stream:
                traceback.print_exc(file=stream)
                error_path = stream.name
            print(f"Comparison failed; detailed diagnostics remain in private log: {error_path}", file=sys.stderr)
        except OSError:
            print("Comparison failed; private diagnostic log unavailable", file=sys.stderr)
        return 1
    print(f"Private exact aggregates: {args.output_root / 'summary.json'}")
    print(f"Disclosure-limited report: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
