"""Synthetic end-to-end comparison; no real clinical fixtures."""
from datetime import datetime, timedelta
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest
import subprocess
import io
from contextlib import redirect_stdout

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from imv_baseline import file_sha256, replay_baseline
from imv_proxy import DEFINITIONS, detect_proxy_episodes, normalize_observations
from compare_imv_definitions import aggregate_proxy_boundaries, run_comparison

ROOT = Path(__file__).resolve().parents[1]
T0 = datetime(2024, 1, 1, 12)


def raw_records(hid, times, devices, *, mode=None):
    return pl.DataFrame({
        "hospitalization_id": [hid] * len(times), "recorded_dttm": [T0 + timedelta(minutes=m) for m in times],
        "device_category": pl.Series(devices, dtype=pl.String), "mode_category": pl.Series(mode or [None] * len(times), dtype=pl.String),
        "tracheostomy": [0] * len(times), "peep_set": [5.0] * len(times), "resp_rate_set": [12.0] * len(times),
        "tidal_volume_set": [400.0] * len(times), "lpm_set": [5.0] + [None] * (len(times)-1),
    })


class ComparisonWorkflow(unittest.TestCase):
    def fixture(self, root):
        source = root / "project"
        (source / "code").mkdir(parents=True)
        (source / "code" / "01_cohort.py").write_bytes((ROOT / "code" / "01_cohort.py").read_bytes())
        wf = raw_records("SYNTHETIC_PRIVATE_H1", [0, 30, 60, 90], ["imv", "imv", "nasal cannula", "room air"])
        hosp = pl.DataFrame({"hospitalization_id": ["SYNTHETIC_PRIVATE_H1"], "discharge_category": ["Home"]})
        cohort = replay_baseline(wf, hosp, source / "code" / "01_cohort.py").with_columns(
            pl.lit("SYNTHETIC_PRIVATE_P1").alias("patient_id"), pl.lit(T0-timedelta(hours=1)).alias("admission_dttm"),
            pl.lit(T0).alias("first_icu_start"), pl.lit(T0+timedelta(hours=4)).alias("first_icu_end"), pl.lit(40).alias("age_at_admission"), pl.lit("Home").alias("discharge_category"),
        )
        cohort.write_parquet(root / "cohort.parquet")
        wf.write_parquet(root / "waterfall.parquet")
        other = raw_records("SYNTHETIC_PRIVATE_H2", [0, 30, 60, 90], [None, None, "Nasal Cannula", "Room Air"],
                            mode=["PRVC", "PRVC", None, None])
        pl.concat([wf, other]).write_parquet(root / "observations.parquet")
        panel = pl.DataFrame({"hospitalization_id": ["SYNTHETIC_PRIVATE_H1", "SYNTHETIC_PRIVATE_H2"],
                              "admission_dttm": [T0-timedelta(hours=1)] * 2, "icu_in_dttm": [T0] * 2,
                              "icu_out_dttm": [T0+timedelta(hours=4)] * 2})
        panel.write_parquet(root / "windows.parquet")
        metadata = {"staged_project": str(source), "cohort": str(root / "cohort.parquet"), "waterfall": str(root / "waterfall.parquet"),
                    "original_observations": str(root / "observations.parquet"), "eligible_windows": str(root / "windows.parquet"),
                    "notebook_sha256": file_sha256(source / "code" / "01_cohort.py"), "original_files_unchanged": True,
                    "parity": {"rows_compared": 1, "mismatches": 0}, "timezone_assumption": "US/Central", "time_semantics": "SYNTHETIC_DATA"}
        path = root / "metadata.json"
        path.write_text(json.dumps(metadata))
        return path

    def test_full_workflow_private_exports_and_redacted_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = self.fixture(root)
            before = file_sha256(root / "cohort.parquet")
            out = root / "private"
            report = root / "report.md"
            console = io.StringIO()
            with redirect_stdout(console):
                summary = run_comparison(metadata, out, report)
            self.assertNotIn("confirmed starts=", console.getvalue())
            self.assertNotIn("confirmed ends=", console.getvalue())
            self.assertNotIn("new candidate admissions=", console.getvalue())
            self.assertEqual(summary["shared"]["B0"]["zero_elapsed_intervals"], 1)
            self.assertEqual(summary["shared"]["A0"]["zero_elapsed_intervals"], 0)
            self.assertEqual(summary["shared"]["C2"]["confirmed_starts"], 1)
            self.assertEqual(summary["discovery"]["C4"]["new_cases"], 1)
            self.assertEqual(summary["discovery"]["C4"]["new_probable_cases"], 1)
            self.assertEqual(before, file_sha256(root / "cohort.parquet"))
            exported = pl.read_parquet(out / "cohorts" / "C2.parquet")
            self.assertEqual(exported["intubation_time"][0], T0)
            self.assertEqual(exported["extubation_time"][0], T0+timedelta(minutes=60))
            self.assertIn("proxy_eligible_end", exported.columns)
            text = report.read_text()
            for sentinel in ("SYNTHETIC_PRIVATE_H1", "SYNTHETIC_PRIVATE_H2", "SYNTHETIC_PRIVATE_P1", "2024-01-01", str(root)):
                self.assertNotIn(sentinel, text)
            if os.name == "posix":
                self.assertEqual(out.stat().st_mode & 0o777, 0o700)
                self.assertEqual((out / "cohorts" / "C2.parquet").stat().st_mode & 0o777, 0o600)

    def test_failed_baseline_validation_stops_before_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = self.fixture(root)
            data = json.loads(path.read_text())
            data["parity"]["mismatches"] = 1
            path.write_text(json.dumps(data))
            with self.assertRaises(RuntimeError) as caught:
                run_comparison(path, root / "private", root / "report.md")
            self.assertNotIn("SYNTHETIC_PRIVATE_H1", str(caught.exception))
            self.assertFalse((root / "report.md").exists())

    def test_gap_segments_do_not_erase_later_observed_cessation(self):
        raw = raw_records("SYNTHETIC_PRIVATE_H1", [0, 30, 240, 270, 300, 330], ["IMV", "IMV", "IMV", "IMV", "Room Air", "Room Air"])
        bounds = pl.DataFrame({"hospitalization_id": ["SYNTHETIC_PRIVATE_H1"], "window_start": [T0-timedelta(hours=1)],
                               "window_end": [T0+timedelta(minutes=400)], "icu_start": [T0]})
        episodes = detect_proxy_episodes(raw, bounds, DEFINITIONS["C2"])
        first = aggregate_proxy_boundaries(episodes).row(0, named=True)
        self.assertEqual(first["start_time"], T0)
        self.assertEqual(first["end_time"], T0+timedelta(minutes=300))
        self.assertTrue(first["interval_has_gap"])
        self.assertEqual(first["segment_count"], 2)

    def test_export_uses_first_cessation_not_final_resumption_end(self):
        raw = raw_records("SYNTHETIC_PRIVATE_H1", [0, 30, 60, 90, 120, 150], ["IMV", "Room Air", "Room Air", "IMV", "Room Air", "Room Air"])
        bounds = pl.DataFrame({"hospitalization_id": ["SYNTHETIC_PRIVATE_H1"], "window_start": [T0-timedelta(hours=1)],
                               "window_end": [T0+timedelta(minutes=200)], "icu_start": [T0]})
        episodes = detect_proxy_episodes(raw, bounds, DEFINITIONS["C2"])
        first = aggregate_proxy_boundaries(episodes).row(0, named=True)
        self.assertEqual(first["end_time"], T0+timedelta(minutes=30))

    def test_all_output_paths_are_checked_before_input_overwrite(self):
        for collision_kind in ("episodes", "provenance"):
            with self.subTest(collision_kind=collision_kind), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                metadata_path = self.fixture(root)
                out = root / "private"
                out.mkdir()
                metadata = json.loads(metadata_path.read_text())
                if collision_kind == "episodes":
                    (out / "episodes").mkdir()
                    protected = out / "episodes" / "C1.parquet"
                    protected.write_bytes((root / "observations.parquet").read_bytes())
                    metadata["original_observations"] = str(protected)
                    metadata_path.write_text(json.dumps(metadata))
                else:
                    protected = out / "provenance.json"
                    protected.write_text(json.dumps(metadata))
                    metadata_path = protected
                before = file_sha256(protected)
                with self.assertRaises(ValueError):
                    run_comparison(metadata_path, out, root / "report.md")
                self.assertEqual(file_sha256(protected), before)

    def test_actual_replay_mismatch_fails_even_if_metadata_claims_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = self.fixture(root)
            cohort = pl.read_parquet(root / "cohort.parquet").with_columns(pl.lit(None, dtype=pl.Datetime("us")).alias("extubation_time"))
            cohort.write_parquet(root / "cohort.parquet")
            with self.assertRaises(RuntimeError) as caught:
                run_comparison(metadata, root / "private", root / "report.md")
            self.assertNotIn("SYNTHETIC_PRIVATE_H1", str(caught.exception))
            self.assertFalse((root / "report.md").exists())

    def test_tracheostomy_after_gap_still_bounds_ett_eligibility(self):
        raw = raw_records("SYNTHETIC_PRIVATE_H1", [0, 30, 240], ["IMV", "IMV", None]).with_columns(pl.Series("tracheostomy", [0, 0, 1]))
        bounds = pl.DataFrame({"hospitalization_id": ["SYNTHETIC_PRIVATE_H1"], "window_start": [T0-timedelta(hours=1)],
                               "window_end": [T0+timedelta(minutes=300)], "icu_start": [T0]})
        episodes = detect_proxy_episodes(raw, bounds, DEFINITIONS["C2"])
        first = aggregate_proxy_boundaries(episodes, normalize_observations(raw)).row(0, named=True)
        self.assertEqual(first["trach_censor_time"], T0+timedelta(minutes=240))
        self.assertIsNone(first["end_time"])

    def test_cohort_nanosecond_keys_join_candidate_microsecond_keys(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = self.fixture(root)
            cohort = pl.read_parquet(root / "cohort.parquet").with_columns(
                *[pl.col(name).cast(pl.Datetime("ns")) for name in ("admission_dttm", "first_icu_start", "first_icu_end")]
            )
            cohort.write_parquet(root / "cohort.parquet")
            summary = run_comparison(metadata, root / "private", root / "report.md", discovery=False)
            self.assertEqual(summary["shared"]["C2"]["confirmed_starts"], 1)
            exported = pl.read_parquet(root / "private" / "cohorts" / "C2.parquet")
            self.assertEqual(exported.schema["first_icu_start"], pl.Datetime("ns"))

    def test_cli_unexpected_schema_error_is_logged_privately(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = self.fixture(root)
            raw = pl.read_parquet(root / "observations.parquet").with_columns(
                pl.Series("device_category", [["SYNTHETIC_PRIVATE_H1"]] * 8, dtype=pl.List(pl.String))
            )
            raw.write_parquet(root / "observations.parquet")
            completed = subprocess.run([sys.executable, str(ROOT / "code" / "compare_imv_definitions.py"),
                                        "--metadata", str(metadata), "--output-root", str(root / "private"), "--report", str(root / "report.md")],
                                       capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            console = completed.stdout + completed.stderr
            self.assertNotIn("Traceback", console)
            self.assertNotIn("SYNTHETIC_PRIVATE_H1", console)
            self.assertTrue(list((root / "private").glob("comparison-error-*.log")))


if __name__ == "__main__":
    unittest.main()
