"""Synthetic records only: proxy timing is not clinical event adjudication."""
from dataclasses import replace
from datetime import datetime, timedelta
from pathlib import Path
import sys
import unittest
from zoneinfo import ZoneInfo

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from imv_proxy import DEFINITIONS, ProxyDefinition, detect_proxy_episodes

T0 = datetime(2024, 1, 1, 12)


def records(entries):
    rows = []
    for entry in entries:
        minute, device, *options = entry
        values = options[0] if options else {}
        row = {"hospitalization_id": "SYNTHETIC_H1", "recorded_dttm": T0 + timedelta(minutes=minute), "device_category": device,
               "tracheostomy": 0, "mode_category": None, "peep_set": None, "resp_rate_set": None,
               "tidal_volume_set": None, "pressure_control_set": None, "lpm_set": None}
        row.update(values)
        rows.append(row)
    return pl.DataFrame(rows, schema_overrides={"mode_category": pl.String, "device_category": pl.String,
        "peep_set": pl.Float64, "resp_rate_set": pl.Float64, "tidal_volume_set": pl.Float64,
        "pressure_control_set": pl.Float64, "lpm_set": pl.Float64})


def windows(end=300):
    return pl.DataFrame({"hospitalization_id": ["SYNTHETIC_H1"], "window_start": [T0 - timedelta(minutes=60)],
                         "window_end": [T0 + timedelta(minutes=end)], "icu_start": [T0]})


def detect(entries, definition="C2", end=300):
    return detect_proxy_episodes(records(entries), windows(end), DEFINITIONS[definition])


class ProxyRules(unittest.TestCase):
    def test_run_span_not_each_adjacent_gap(self):
        result = detect([(-60, "Room Air"), (0, "IMV"), (5, "IMV"), (15, "IMV"), (30, "IMV"),
                         (60, "Nasal Cannula"), (75, "Room Air"), (90, "Room Air")], "C1")
        row = result.row(0, named=True)
        self.assertEqual(row["start_time"], T0)
        self.assertEqual(row["start_confirmation_time"], T0 + timedelta(minutes=30))
        self.assertEqual(row["end_time"], T0 + timedelta(minutes=60))
        self.assertEqual(row["end_confirmation_time"], T0 + timedelta(minutes=90))
        self.assertFalse(row["left_censored"])

    def test_immediate_start_retains_a_brief_explicit_episode(self):
        rows = [(0, "IMV"), (5, "Nasal Cannula"), (35, "Room Air")]
        balanced = detect(rows, "C1").row(0, named=True)
        sensitive = detect(rows, "C2").row(0, named=True)
        self.assertFalse(balanced["start_confirmed"])
        self.assertIsNone(balanced["end_time"])
        self.assertTrue(sensitive["start_confirmed"])
        self.assertEqual(sensitive["end_time"], T0 + timedelta(minutes=5))

    def test_controlled_to_pressure_support_mode_is_not_an_end(self):
        result = detect([(0, "IMV", {"mode_category": "PRVC"}),
                         (30, "IMV", {"mode_category": "Pressure Support/CPAP"}), (60, "IMV")])
        self.assertIsNone(result["end_time"][0])
        self.assertEqual(result["censor_reason"][0], "window_end")

    def test_niv_settings_never_override_explicit_device(self):
        setting = {"mode_category": "Pressure Control", "peep_set": 5, "resp_rate_set": 12, "tidal_volume_set": 400}
        result = detect([(0, "NIPPV", setting), (30, "NIPPV", setting)], "C4")
        self.assertEqual(result.height, 0)

    def test_oxygen_flow_on_imv_is_not_an_end(self):
        result = detect([(0, "IMV", {"lpm_set": 5}), (30, "IMV", {"lpm_set": 5})])
        self.assertTrue(result["start_confirmed"][0])
        self.assertIsNone(result["end_time"][0])

    def test_missing_device_can_be_probable_not_direct(self):
        setting = {"mode_category": "Pressure-Regulated Volume Control", "peep_set": 5, "resp_rate_set": 12, "tidal_volume_set": 400}
        rows = [(0, None, setting), (30, None, setting), (60, "Nasal Cannula"), (90, "Room Air")]
        self.assertEqual(detect(rows, "C2").height, 0)
        row = detect(rows, "C4").row(0, named=True)
        self.assertTrue(row["start_confirmed"])
        self.assertEqual(row["start_evidence"], "probable")
        self.assertEqual(row["start_confirmation_time"], T0 + timedelta(minutes=30))

    def test_peep_alone_or_nonexclusive_mode_does_not_infer_imv(self):
        peep = {"peep_set": 5}
        ps = {"mode_category": "Pressure Support/CPAP", "peep_set": 5, "resp_rate_set": 12, "tidal_volume_set": 400}
        self.assertEqual(detect([(0, None, peep), (30, None, peep)], "C4").height, 0)
        self.assertEqual(detect([(0, None, ps), (30, None, ps)], "C4").height, 0)

    def test_unknown_rows_do_not_confirm_or_end(self):
        result = detect([(0, "IMV"), (15, None), (30, "IMV"), (60, None)], "C1")
        self.assertTrue(result["start_confirmed"][0])
        self.assertEqual(result["start_confirmation_time"][0], T0 + timedelta(minutes=30))
        self.assertIsNone(result["end_time"][0])
        self.assertTrue(result["had_unknown_observations"][0])

    def test_documentation_gap_censors_without_inventing_an_end(self):
        result = detect([(0, "IMV"), (30, "IMV"), (240, "IMV")], "C1", end=300)
        self.assertEqual(result.height, 2)
        self.assertEqual(result["censor_reason"][0], "documentation_gap")
        self.assertIsNone(result["end_time"][0])
        self.assertFalse(result["start_confirmed"][1])

    def test_duplicate_timestamp_is_not_two_observations(self):
        result = detect([(0, "IMV"), (0, "IMV"), (30, "Room Air"), (60, "Room Air")], "C1")
        self.assertFalse(result["start_confirmed"][0])

    def test_scaffold_does_not_supply_confirmation(self):
        result = detect([(0, "IMV", {"is_scaffold": False}), (30, "IMV", {"is_scaffold": True}),
                         (60, "Room Air", {"is_scaffold": False})], "C1")
        self.assertFalse(result["start_confirmed"][0])

    def test_opposing_same_time_devices_are_unresolved(self):
        result = detect([(0, "IMV"), (0, "NIPPV"), (30, "Room Air")], "C4")
        self.assertEqual(result.height, 0)

    def test_new_tracheostomy_is_censoring_not_extubation(self):
        row = detect([(0, "IMV"), (30, "IMV", {"tracheostomy": 1})]).row(0, named=True)
        self.assertTrue(row["start_confirmed"])
        self.assertIsNone(row["end_time"])
        self.assertEqual(row["censor_reason"], "new_tracheostomy")
        self.assertEqual(detect([(0, "IMV", {"tracheostomy": 1})]).height, 0)

    def test_unknown_tracheostomy_is_not_assumed_absent(self):
        self.assertEqual(detect([(0, "IMV", {"tracheostomy": None}), (30, "IMV", {"tracheostomy": None})]).height, 0)

    def test_resumption_changes_30_and_60_minute_end_rules(self):
        rows = [(0, "IMV"), (30, "Nasal Cannula"), (60, "Room Air"), (90, "IMV"),
                (120, "Nasal Cannula"), (150, "Room Air")]
        sensitive = detect(rows, "C2", end=180)
        conservative = detect(rows, "C3", end=180)
        self.assertEqual(sensitive.height, 2)
        self.assertEqual(sensitive["end_time"][0], T0 + timedelta(minutes=30))
        self.assertEqual(conservative.height, 1)
        self.assertIsNone(conservative["end_time"][0])

    def test_event_and_sat_sbt_labels_are_not_inputs(self):
        raw = records([(0, "IMV"), (30, "IMV"), (60, "Room Air"), (90, "Room Air")])
        tagged = raw.with_columns(pl.lit(False).alias("is_intubated"), pl.lit(True).alias("extubation_happened"),
                                 pl.lit(1).alias("consensus_ett_extubation"), pl.lit(1).alias("sbt_delivery_performed"))
        left = detect_proxy_episodes(raw, windows(), DEFINITIONS["C2"])
        right = detect_proxy_episodes(tagged, windows(), DEFINITIONS["C2"])
        self.assertTrue(left.equals(right))

    def test_clif3_device_aliases(self):
        row = detect([(0, "imv"), (30, "imv"), (60, "nasal_cannula"), (90, "room_air")]).row(0, named=True)
        self.assertEqual(row["end_time"], T0 + timedelta(minutes=60))

    def test_empty_inputs_and_bad_schema(self):
        empty = records([(0, "IMV")]).head(0)
        self.assertEqual(detect_proxy_episodes(empty, windows(), DEFINITIONS["C2"]).height, 0)
        with self.assertRaises(ValueError):
            detect_proxy_episodes(empty.drop("tracheostomy"), windows(), DEFINITIONS["C2"])
        with self.assertRaises(ValueError):
            detect_proxy_episodes(records([(0, "IMV")]).with_columns(pl.col("recorded_dttm").cast(pl.String)), windows(), DEFINITIONS["C2"])
        with self.assertRaises(ValueError):
            replace(DEFINITIONS["C1"], start_minutes=-1)

    def test_aware_time_uses_elapsed_utc_through_dst_fold(self):
        zone = ZoneInfo("America/Chicago")
        first = datetime(2024, 11, 3, 1, 0, tzinfo=zone, fold=0)
        second = datetime(2024, 11, 3, 1, 30, tzinfo=zone, fold=1)
        raw = records([(0, "IMV"), (30, "IMV")]).with_columns(pl.Series("recorded_dttm", [first, second]))
        bounds = pl.DataFrame({"hospitalization_id": ["SYNTHETIC_H1"], "window_start": [first - timedelta(hours=2)],
                               "window_end": [second + timedelta(hours=2)], "icu_start": [first]})
        result = detect_proxy_episodes(raw, bounds, replace(DEFINITIONS["C1"], start_minutes=60))
        row = result.row(0, named=True)
        self.assertTrue(row["start_confirmed"])
        self.assertEqual((row["start_confirmation_time"] - row["start_time"]).total_seconds(), 90 * 60)


if __name__ == "__main__":
    unittest.main()
