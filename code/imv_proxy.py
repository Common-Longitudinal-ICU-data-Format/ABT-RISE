"""Independent EHR proxies for ETT-IMV onset/cessation, not proven airway events.

Explicit device states take precedence. Ventilator-setting inference is separately
labeled probable. Dedicated airway-event and SAT/SBT labels are never inputs.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import math
from typing import Any

import polars as pl


@dataclass(frozen=True)
class ProxyDefinition:
    name: str = "custom"
    start_minutes: float = 30
    end_minutes: float = 30
    start_observations: int = 2
    end_observations: int = 2
    evidence_assisted: bool = False
    inferred_start_minutes: float = 30
    max_gap_minutes: float = 120

    def __post_init__(self) -> None:
        for value in (self.start_minutes, self.end_minutes, self.inferred_start_minutes):
            if not math.isfinite(value) or value < 0:
                raise ValueError("Confirmation periods must be finite and nonnegative")
        if not math.isfinite(self.max_gap_minutes) or self.max_gap_minutes <= 0:
            raise ValueError("Maximum evidence gap must be finite and positive")
        if self.start_observations < 1 or self.end_observations < 2:
            raise ValueError("Start requires observations; end requires at least two")


DEFINITIONS = {
    "C1": ProxyDefinition(name="C1"),
    "C2": ProxyDefinition(name="C2", start_minutes=0, start_observations=1),
    "C3": ProxyDefinition(name="C3", start_minutes=0, start_observations=1, end_minutes=60),
    "C4": ProxyDefinition(name="C4", start_minutes=0, start_observations=1, evidence_assisted=True),
    "C0": ProxyDefinition(name="C0", start_minutes=0, end_minutes=0),
}

OFF_DEVICES = {"nippv", "cpap", "high flow nc", "face mask", "trach collar", "nasal cannula", "room air", "t piece", "blow by"}
CONTROLLED_MODES = {"assist control volume control", "pressure control", "pressure regulated volume control", "simv"}
DEVICE_ALIASES = {"hfnc": "high flow nc"}
MODE_ALIASES = {"acvc": "assist control volume control", "prvc": "pressure regulated volume control", "ps or cpap": "pressure support/cpap"}
TIME_FIELDS = {
    "window_start", "window_end", "icu_start", "start_time", "start_confirmation_time", "start_lower_bound",
    "end_time", "end_confirmation_time", "end_lower_bound", "end_upper_bound", "last_imv_time", "censor_time",
}


def event_schema(timezone: str | None = None) -> dict[str, pl.DataType]:
    schema: dict[str, pl.DataType] = {
        "hospitalization_id": pl.String, "definition": pl.String, "episode_number": pl.Int64,
        "window_start": pl.Datetime("us", timezone), "window_end": pl.Datetime("us", timezone), "icu_start": pl.Datetime("us", timezone),
        "start_time": pl.Datetime("us", timezone), "start_confirmation_time": pl.Datetime("us", timezone),
        "start_confirmed": pl.Boolean, "start_evidence": pl.String, "start_lower_bound": pl.Datetime("us", timezone), "left_censored": pl.Boolean,
        "end_time": pl.Datetime("us", timezone), "end_confirmation_time": pl.Datetime("us", timezone), "end_confirmed": pl.Boolean,
        "end_lower_bound": pl.Datetime("us", timezone), "end_upper_bound": pl.Datetime("us", timezone), "last_imv_time": pl.Datetime("us", timezone),
        "censor_time": pl.Datetime("us", timezone), "censor_reason": pl.String,
        "had_unknown_observations": pl.Boolean, "had_conflicting_observations": pl.Boolean,
        "used_probable_evidence": pl.Boolean, "has_icu_imv_evidence": pl.Boolean,
        "start_support_points": pl.Int64, "end_support_points": pl.Int64, "duration_minutes": pl.Float64,
    }
    return schema


def _text(column: str, aliases: dict[str, str]) -> pl.Expr:
    expression = pl.col(column).cast(pl.String).str.strip_chars().str.to_lowercase().str.replace_all("[_-]", " ")
    return expression.replace(aliases)


def normalize_observations(observations: pl.DataFrame) -> pl.DataFrame:
    """Normalize vocabulary, classify contemporaneous evidence, collapse timestamps.

The caller must supply original observations. Excluding scaffolds cannot reverse
imputation on original rows in an already-waterfalled export.
"""
    required = {"hospitalization_id", "recorded_dttm", "device_category", "tracheostomy"}
    if not required.issubset(observations.columns):
        raise ValueError("Respiratory observations are missing required columns")
    dtype = observations.schema["recorded_dttm"]
    if not isinstance(dtype, pl.Datetime):
        raise ValueError("Respiratory observation times must be datetime values")
    frame = observations.filter(pl.col("hospitalization_id").is_not_null() & pl.col("recorded_dttm").is_not_null())
    if "is_scaffold" in frame.columns:
        frame = frame.filter(~pl.col("is_scaffold").fill_null(False))
    if "mode_category" not in frame.columns:
        frame = frame.with_columns(pl.lit(None, dtype=pl.String).alias("mode_category"))
    numbers = ("peep_set", "resp_rate_set", "tidal_volume_set", "pressure_control_set")
    for name in numbers:
        if name not in frame.columns:
            frame = frame.with_columns(pl.lit(None, dtype=pl.Float64).alias(name))
    timestamp = pl.col("recorded_dttm")
    if dtype.time_zone is not None:
        timestamp = timestamp.dt.convert_time_zone("UTC")
    frame = frame.with_columns(
        pl.col("hospitalization_id").cast(pl.String), timestamp.cast(pl.Datetime("us", "UTC" if dtype.time_zone else None)),
        _text("device_category", DEVICE_ALIASES).alias("_device"), _text("mode_category", MODE_ALIASES).alias("_mode"),
        pl.col("tracheostomy").cast(pl.Float64, strict=False).alias("_trach"),
        *[pl.col(name).cast(pl.Float64, strict=False).fill_nan(None) for name in numbers],
    ).with_columns(pl.when(pl.col("_trach").is_in([0.0, 1.0])).then(pl.col("_trach")).otherwise(None).alias("_trach"))
    setting_evidence = (
        pl.col("peep_set").is_finite() & (pl.col("peep_set") >= 0)
        & pl.col("resp_rate_set").is_finite() & (pl.col("resp_rate_set") > 0)
        & ((pl.col("tidal_volume_set").is_finite() & (pl.col("tidal_volume_set") > 0))
           | ((pl.col("_mode") == "pressure control") & pl.col("pressure_control_set").is_finite() & (pl.col("pressure_control_set") > 0)))
    )
    frame = frame.with_columns(
        pl.when(pl.col("_device") == "imv").then(1)
        .when(pl.col("_device").is_in(OFF_DEVICES)).then(0)
        .when(pl.col("_device").is_null() & pl.col("_mode").is_in(CONTROLLED_MODES) & setting_evidence).then(2)
        .otherwise(-1).alias("_signal")
    )
    grouped = frame.group_by("hospitalization_id", "recorded_dttm").agg(
        pl.col("_signal").unique().alias("_signals"), pl.col("_trach").max(), pl.len().alias("_original_rows"),
    )
    on = pl.col("_signals").list.contains(1)
    off = pl.col("_signals").list.contains(0)
    probable = pl.col("_signals").list.contains(2)
    return grouped.with_columns(
        (on & off).alias("_conflict"),
        pl.when(on & off).then(-1).when(on).then(1).when(off).then(0).when(probable).then(2).otherwise(-1).alias("_signal"),
    ).select("hospitalization_id", "recorded_dttm", "_signal", "_trach", "_conflict", "_original_rows").sort("hospitalization_id", "recorded_dttm")


def _normalize_windows(windows: pl.DataFrame, observation_timezone: str | None) -> pl.DataFrame:
    required = {"hospitalization_id", "window_start", "window_end"}
    if not required.issubset(windows.columns):
        raise ValueError("Eligible windows are missing required columns")
    if "icu_start" not in windows.columns:
        windows = windows.with_columns(pl.col("window_start").alias("icu_start"))
    columns = ("window_start", "window_end", "icu_start")
    expressions = [pl.col("hospitalization_id").cast(pl.String)]
    for name in columns:
        dtype = windows.schema[name]
        if not isinstance(dtype, pl.Datetime) or (dtype.time_zone is None) != (observation_timezone is None):
            raise ValueError("Observation and window datetime awareness must agree")
        value = pl.col(name)
        if dtype.time_zone is not None:
            value = value.dt.convert_time_zone("UTC")
        expressions.append(value.cast(pl.Datetime("us", observation_timezone)))
    windows = windows.with_columns(expressions)
    if windows.filter(pl.any_horizontal([pl.col(name).is_null() for name in (*columns, "hospitalization_id")]) | (pl.col("window_start") > pl.col("window_end"))).height:
        raise ValueError("Eligible windows have missing or reversed bounds")
    return windows.unique(subset=["hospitalization_id", "window_start", "window_end", "icu_start"]).sort("hospitalization_id", "icu_start", "window_end")


def _minutes(later: datetime, earlier: datetime) -> float:
    return (later - earlier).total_seconds() / 60


def detect_proxy_episodes(
    observations: pl.DataFrame, windows: pl.DataFrame, definition: ProxyDefinition = DEFINITIONS["C1"],
) -> pl.DataFrame:
    """Return confirmed and unconfirmed ETT-IMV segments under one declared policy."""
    normalized = normalize_observations(observations)
    return detect_normalized_episodes(normalized, windows, definition)


def detect_normalized_episodes(normalized: pl.DataFrame, windows: pl.DataFrame, definition: ProxyDefinition) -> pl.DataFrame:
    """Execute on a reusable classified table produced by normalize_observations."""
    timezone = normalized.schema["recorded_dttm"].time_zone
    windows = _normalize_windows(windows, timezone)
    schema = event_schema(timezone)
    groups = normalized.partition_by("hospitalization_id", as_dict=True, maintain_order=True)
    output: list[dict[str, Any]] = []
    for window in windows.iter_rows(named=True):
        group = groups.get((window["hospitalization_id"],))
        if group is None:
            continue
        rows = group.filter(pl.col("recorded_dttm").is_between(window["window_start"], window["window_end"]))
        pending: dict[str, Any] | None = None
        off_start: datetime | None = None
        off_points = 0
        last_known: datetime | None = None
        previous_off: datetime | None = None
        last_trach: float | None = None
        episode_number = 0
        trach_seen = False

        def emit(reason: str | None, censor_time: datetime | None = None) -> None:
            nonlocal pending, off_start, off_points
            if pending is None:
                return
            pending["end_support_points"] = off_points
            pending["end_upper_bound"] = off_start
            pending["end_lower_bound"] = pending["last_imv_time"]
            pending["censor_reason"] = reason
            pending["censor_time"] = censor_time
            if pending["end_time"] is not None and pending["start_confirmed"]:
                pending["duration_minutes"] = _minutes(pending["end_time"], pending["start_time"])
            output.append({name: pending.get(name) for name in schema})
            pending = None
            off_start = None
            off_points = 0

        for row in rows.iter_rows(named=True):
            time = row["recorded_dttm"]
            signal = row["_signal"]
            if signal == 2 and not definition.evidence_assisted:
                signal = -1
            if last_known is not None and _minutes(time, last_known) > definition.max_gap_minutes:
                emit("documentation_gap", last_known)
                last_known = None
                last_trach = None
            if row["_trach"] is not None:
                last_trach = row["_trach"]
            if last_trach == 1:
                emit("new_tracheostomy", time)
                trach_seen = True
            if trach_seen:
                continue
            if signal < 0:
                if pending is not None:
                    pending["had_unknown_observations"] = True
                    pending["had_conflicting_observations"] |= row["_conflict"]
                continue
            if signal in (1, 2):
                if last_trach != 0:
                    # Unknown route is not silently treated as endotracheal IMV.
                    if pending is not None:
                        pending["had_unknown_observations"] = True
                    continue
                last_known = time
                if pending is not None and not pending["start_confirmed"] and pending["start_evidence"] == "probable" and signal == 1:
                    # An unconfirmed inferred onset does not delay explicit IMV entry.
                    pending["start_time"] = time
                    pending["start_evidence"] = "direct"
                    pending["start_support_points"] = 0
                if pending is None:
                    episode_number += 1
                    pending = {
                        "hospitalization_id": window["hospitalization_id"], "definition": definition.name,
                        "episode_number": episode_number, "window_start": window["window_start"], "window_end": window["window_end"], "icu_start": window["icu_start"],
                        "start_time": time, "start_confirmation_time": None, "start_confirmed": False,
                        "start_evidence": "direct" if signal == 1 else "probable", "start_lower_bound": previous_off,
                        "left_censored": previous_off is None, "end_time": None, "end_confirmation_time": None,
                        "end_confirmed": False, "last_imv_time": time, "censor_reason": None,
                        "had_unknown_observations": False, "had_conflicting_observations": False,
                        "used_probable_evidence": signal == 2, "has_icu_imv_evidence": False,
                        "start_support_points": 0, "end_support_points": 0, "duration_minutes": None,
                    }
                off_start = None
                off_points = 0
                pending["last_imv_time"] = time
                pending["used_probable_evidence"] |= signal == 2
                pending["has_icu_imv_evidence"] |= time >= window["icu_start"]
                if not pending["start_confirmed"]:
                    pending["start_support_points"] += 1
                    inferred = pending["start_evidence"] == "probable"
                    span = definition.inferred_start_minutes if inferred else definition.start_minutes
                    count = 2 if inferred else definition.start_observations
                    if pending["start_support_points"] >= count and _minutes(time, pending["start_time"]) >= span:
                        pending["start_confirmed"] = True
                        pending["start_confirmation_time"] = time
            else:
                last_known = time
                previous_off = time
                if pending is None:
                    continue
                if not pending["start_confirmed"]:
                    off_start = time
                    off_points = 1
                    emit("unconfirmed_start", time)
                    continue
                if off_start is None:
                    off_start = time
                off_points += 1
                if off_points >= definition.end_observations and _minutes(time, off_start) >= definition.end_minutes:
                    pending["end_time"] = off_start
                    pending["end_confirmed"] = True
                    pending["end_confirmation_time"] = time
                    emit(None)
        emit("window_end", window["window_end"])
    return pl.DataFrame(output, schema=schema, strict=False) if output else pl.DataFrame(schema=schema)
