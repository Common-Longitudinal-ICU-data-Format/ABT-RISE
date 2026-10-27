"""Synthetic characterization of the unmodified cohort event cell."""
from datetime import datetime, timedelta
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from imv_baseline import (
    assert_baseline_parity,
    replay_baseline,
    stage_baseline,
    file_sha256,
    run_baseline,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "code" / "01_cohort.py"
T0 = datetime(2024, 1, 1, 12)


def observations(devices, flow=None):
    return pl.DataFrame({
        "hospitalization_id": ["SYNTHETIC_H1"] * len(devices),
        "recorded_dttm": [T0 + timedelta(minutes=i * 60) for i in range(len(devices))],
        "device_category": devices,
        "lpm_set": pl.Series(flow or [None] * len(devices), dtype=pl.Float64),
    })


def hospitals():
    return pl.DataFrame({"hospitalization_id": ["SYNTHETIC_H1"], "discharge_category": ["Home"]})


class BaselineCharacterization(unittest.TestCase):
    def test_flow_shortcut_preserves_zero_time_end(self):
        result = replay_baseline(observations(["imv"], [5.0]), hospitals(), SOURCE)
        row = result.row(0, named=True)
        self.assertEqual(row["intubation_time"], T0)
        self.assertEqual(row["extubation_time"], T0)
        self.assertEqual(row["extubation_method"], "2-row IMV+LPM")

    def test_device_end_is_last_imv_not_first_off(self):
        result = replay_baseline(observations(["imv", "imv", "nasal cannula", "room air"]), hospitals(), SOURCE)
        self.assertEqual(result["extubation_time"][0], T0 + timedelta(hours=1))
        self.assertEqual(result["extubation_method"][0], "2-row IMV")

    def test_single_off_reading_is_unconfirmed(self):
        result = replay_baseline(observations(["imv", "nasal cannula"]), hospitals(), SOURCE)
        self.assertIsNone(result["extubation_time"][0])
        self.assertIn("Unconfirmed", result["no_extubation_reason"][0])

    def test_missing_cell_and_wrong_source_hash_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = Path(tmp) / "notebook.py"
            missing.write_text("x = 1\n")
            with self.assertRaisesRegex(ValueError, "event cell"):
                replay_baseline(observations(["imv"]), hospitals(), missing)
        with self.assertRaisesRegex(ValueError, "source hash"):
            replay_baseline(observations(["imv"]), hospitals(), SOURCE, expected_sha256="wrong")

    def test_parity_detects_null_and_nonnull_mismatch_without_rows(self):
        result = replay_baseline(observations(["imv"], [5.0]), hospitals(), SOURCE)
        self.assertEqual(assert_baseline_parity(result, result)["rows_compared"], 1)
        changed = result.with_columns(pl.lit(None, dtype=pl.Datetime("us")).alias("extubation_time"))
        with self.assertRaises(RuntimeError) as caught:
            assert_baseline_parity(result, changed)
        self.assertNotIn("SYNTHETIC_H1", str(caught.exception))
        self.assertNotIn("2024-01", str(caught.exception))

    def test_stage_is_byte_identical_and_does_not_touch_original(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "source"
            (root / "code").mkdir(parents=True)
            (root / "code" / "01_cohort.py").write_bytes(SOURCE.read_bytes())
            (root / "code" / "helper.py").write_text("# synthetic helper\n")
            (root / "clif_config.json").write_text(json.dumps({"data_directory": "data", "timezone": "US/Central"}))
            (root / "output_phi").mkdir()
            original = root / "output_phi" / "cohort.parquet"
            original.write_bytes(b"SYNTHETIC ORIGINAL")
            before = original.read_bytes()
            staged = stage_baseline(root, root / "output_phi" / "imv_comparison")
            self.assertEqual(file_sha256(root / "code" / "01_cohort.py"), file_sha256(staged / "code" / "01_cohort.py"))
            self.assertEqual(original.read_bytes(), before)
            self.assertNotEqual(staged, root)
            renamed_cache = Path(tmp) / "archived-waterfall.parquet"
            renamed_cache.write_bytes(b"SYNTHETIC CACHE BYTES")
            with_cache = stage_baseline(root, root / "output_phi" / "imv_comparison", waterfall_cache=renamed_cache)
            self.assertEqual((with_cache / "output_phi" / "resp_waterfall_cohort.parquet").read_bytes(), renamed_cache.read_bytes())
            cfg = json.loads((staged / "clif_config.json").read_text())
            self.assertEqual(Path(cfg["data_directory"]), (root / "data").resolve())
            if os.name == "posix":
                self.assertEqual(staged.stat().st_mode & 0o777, 0o700)
                self.assertEqual((staged / "clif_config.json").stat().st_mode & 0o777, 0o600)

    def test_method_label_difference_requires_proven_equal_branch_times(self):
        frame = observations(["imv", "imv", "nasal cannula", "room air"], [None, 5.0, None, None])
        result, ties = replay_baseline(frame, hospitals(), SOURCE, return_branch_candidates=True)
        alternate = result.with_columns(
            pl.when(pl.col("extubation_method") == "2-row IMV")
            .then(pl.lit("2-row IMV+LPM")).otherwise(pl.lit("2-row IMV")).alias("extubation_method")
        )
        with self.assertRaises(RuntimeError):
            assert_baseline_parity(result, alternate)
        parity = assert_baseline_parity(result, alternate, method_ties=ties)
        self.assertEqual(parity["method_tie_disagreements"], 1)
        changed_time = alternate.with_columns((pl.col("extubation_time") + pl.duration(minutes=1)).alias("extubation_time"))
        with self.assertRaises(RuntimeError):
            assert_baseline_parity(result, changed_time, method_ties=ties)

    def test_method_difference_without_two_equal_candidates_still_fails(self):
        result, ties = replay_baseline(observations(["imv"], [5.0]), hospitals(), SOURCE, return_branch_candidates=True)
        alternate = result.with_columns(pl.lit("2-row IMV").alias("extubation_method"))
        self.assertEqual(ties.height, 0)
        with self.assertRaises(RuntimeError):
            assert_baseline_parity(result, alternate, method_ties=ties)

    def test_notebook_default_outputs_use_private_working_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve() / "source"
            (root / "code").mkdir(parents=True)
            (root / "data").mkdir()
            (root / "data" / "clif_respiratory_support.parquet").write_bytes(b"SYNTHETIC_INPUT")
            (root / "code" / "01_cohort.py").write_bytes(SOURCE.read_bytes())
            (root / "code" / "helper.py").write_text("# synthetic helper\n")
            (root / "clif_config.json").write_text(json.dumps({"data_directory": "data", "timezone": "US/Central"}))
            output = root / "output_phi" / "imv_comparison"
            original_cwd = Path.cwd()

            def fake_run():
                staged = Path.cwd()
                self.assertTrue(staged.is_relative_to(output), "Notebook working directory must be private staging")
                default_output = staged / "output"
                default_output.mkdir()
                (default_output / "synthetic.txt").write_text("SYNTHETIC_DEFAULT_OUTPUT")
                resp = observations(["imv", "imv"], [5.0, None]).with_columns(pl.lit(0).alias("tracheostomy"))
                hosp = hospitals().with_columns(
                    pl.lit("SYNTHETIC_P1").alias("patient_id"),
                    pl.lit(T0 - timedelta(hours=1)).alias("admission_dttm"),
                    pl.lit(T0 + timedelta(days=2)).alias("discharge_dttm"),
                    pl.lit(40).alias("age_at_admission"),
                )
                cohort = replay_baseline(resp, hosp, SOURCE)
                (staged / "output_phi").mkdir()
                cohort.write_parquet(staged / "output_phi" / "cohort.parquet")
                resp.write_parquet(staged / "output_phi" / "resp_waterfall_cohort.parquet")
                icu = pl.DataFrame({"hospitalization_id": ["SYNTHETIC_H1"], "icu_in_dttm": [T0], "icu_out_dttm": [T0 + timedelta(days=1)]})
                return [], {"cohort_df": cohort, "resp_wf": resp, "hosp_final": hosp, "hosp_24h": hosp, "rs_all_df": resp, "icu_24h": icu}

            spec = SimpleNamespace(loader=SimpleNamespace(exec_module=lambda module: setattr(module, "app", SimpleNamespace(run=fake_run))))
            with patch("imv_baseline.importlib.util.spec_from_file_location", return_value=spec), patch("imv_baseline.importlib.util.module_from_spec", return_value=SimpleNamespace()):
                metadata = run_baseline(root, output)
            self.assertTrue(metadata["original_files_unchanged"])
            self.assertEqual(Path.cwd(), original_cwd)
            self.assertFalse((root / "output").exists())
            self.assertTrue((Path(metadata["staged_project"]) / "output" / "synthetic.txt").exists())


if __name__ == "__main__":
    unittest.main()
