import marimo

__generated_with = "0.20.2"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    import polars as pl
    import pandas as pd
    import numpy as np
    import json
    from pathlib import Path
    from tqdm import tqdm
    from helper import plot_consort, plot_upset

    # Read site-specific config (site name, paths, etc.)
    _config_path = Path(__file__).parent.parent / "clif_config.json"
    with open(_config_path) as _f:
        _cfg = json.load(_f)

    SITE_NAME = _cfg["site_name"]
    DATA_DIR = _cfg["data_directory"]
    TIMEZONE = _cfg.get("timezone", None)
    OUTPUT_PHI = Path(__file__).parent.parent / "output_phi"
    OUTPUT_SBT = OUTPUT_PHI / "sbt_both_stabilities"        # row-level PHI data
    OUTPUT_SBT.mkdir(parents=True, exist_ok=True)
    OUTPUT_SHARE = Path(__file__).parent.parent / "output_to_share" / "sbt_both_stabilities"
    OUTPUT_SHARE.mkdir(parents=True, exist_ok=True)

    print(f"Site: {SITE_NAME}")
    print(f"PHI output:   {OUTPUT_SBT}")
    print(f"Share output: {OUTPUT_SHARE}")
    return (
        DATA_DIR,
        OUTPUT_PHI,
        OUTPUT_SBT,
        OUTPUT_SHARE,
        Path,
        SITE_NAME,
        TIMEZONE,
        json,
        np,
        pd,
        pl,
        plot_consort,
        plot_upset,
        tqdm,
    )


@app.cell
def _(OUTPUT_PHI, pl):
    wide = pl.read_parquet(OUTPUT_PHI / "wide_dataset.parquet")
    print(f"Wide dataset: {wide.height:,} rows x {wide.width} cols")

    cohort = pl.read_parquet(OUTPUT_PHI / "cohort.parquet")
    print(f"Cohort: {cohort.height:,} hospitalizations")
    return cohort, wide


@app.cell
def _(pl, wide):
    df = wide.with_columns(
        pl.col("recorded_dttm").cast(pl.Datetime("us")),
        pl.col("first_icu_start").cast(pl.Datetime("us")),
        pl.col("first_icu_end").cast(pl.Datetime("us")),
    ).sort("hospitalization_id", "recorded_dttm")
    print(f"Preprocessed: {df.height:,} rows, {df['hosp_id_icu_day'].n_unique():,} hosp-day keys")
    df = df.with_columns(
        pl.col("nee").fill_null(0.0).alias("nee_filled"),
    ).with_columns(
        (pl.col("nee_filled") <= 0.2).cast(pl.Int32).alias("hemodynamic_stable"),
    ).drop("nee_filled")
    df = df.with_columns(
        (
            (pl.col("fio2_set").fill_null(0.21) <= 0.5)
            & (pl.col("peep_set").fill_null(0.0) <= 8)
            & (pl.col("spo2").fill_null(100.0) >= 88)
        ).cast(pl.Int32).alias("respiratory_stable"),
    )
    print(f"Hemodynamic stability: {df.filter(pl.col('hemodynamic_stable') == 1).height:,} / {df.height:,} rows stable")
    print(f"Respiratory stability: {df.filter(pl.col('respiratory_stable') == 1).height:,} / {df.height:,} rows stable")
    _both = df.filter((pl.col('hemodynamic_stable') == 1) & (pl.col('respiratory_stable') == 1)).height
    print(f"Both stable: {_both:,} / {df.height:,} rows")
    return (df,)


@app.cell
def _(df, pl):
    # All unique patient-days (needed so days with zero overnight rows still appear)
    _all_keys = df.select("hosp_id_icu_day").unique()

    # Filter to overnight window (22:00–06:00), same as the 6h eligibility window
    _overnight = df.filter(
        (pl.col("recorded_dttm") >= pl.col("icu_day_date") - pl.duration(hours=2))
        & (pl.col("recorded_dttm") <= pl.col("icu_day_date") + pl.duration(hours=6))
    )

    # Per day: check each condition within the overnight window only
    # SBT eligibility: IMV + no paralytics (no sedation requirement)
    _overnight_flags = (
        _overnight.group_by("hosp_id_icu_day")
        .agg(
            (pl.col("max_paralytics") > 0).any().alias("has_paralytics"),
            (pl.col("device_category").str.to_lowercase() == "imv").any().alias("has_imv"),
            (pl.col("respiratory_stable") == 1).any().alias("has_resp_stable"),
            (pl.col("hemodynamic_stable") == 1).any().alias("has_hemo_stable"),
        )
    )

    # Left-join onto all days; days with no overnight data → null before fill
    _day_flags = _all_keys.join(_overnight_flags, on="hosp_id_icu_day", how="left")

    # Flag days that had ANY rows in the overnight window (non-null before fill)
    _day_flags = _day_flags.with_columns(
        pl.col("has_paralytics").is_not_null().alias("has_overnight_data"),
    )

    _day_flags = _day_flags.with_columns(
        pl.col("has_paralytics").fill_null(False),
        pl.col("has_imv").fill_null(False),
        pl.col("has_resp_stable").fill_null(False),
        pl.col("has_hemo_stable").fill_null(False),
    )

    # Hierarchical failure reason: no overnight data 1st, then no IMV, then paralytics, then <6h
    _day_flags = _day_flags.with_columns(
        pl.when(~pl.col("has_overnight_data")).then(pl.lit("No overnight data"))
        .when(~pl.col("has_imv")).then(pl.lit("No IMV"))
        .when(pl.col("has_paralytics")).then(pl.lit("Paralytics present"))
        .otherwise(pl.lit("< 6h continuous window"))
        .alias("eligibility_failure_reason")
    )

    failure_reason_df = _day_flags.select("hosp_id_icu_day", "eligibility_failure_reason")

    # Independent boolean failure flags (one row per IMV day, multiple True possible)
    failure_flags_df = (
        _day_flags.filter(pl.col("has_overnight_data") & pl.col("has_imv"))
        .select(
            "hosp_id_icu_day",
            pl.col("has_paralytics"),
            (~pl.col("has_resp_stable")).alias("resp_unstable"),
            (~pl.col("has_hemo_stable")).alias("hemo_unstable"),
        )
    )

    # CONSORT step counts: total → overnight → IMV (remaining go to eligibility check)
    _total = _day_flags.height
    _n_overnight = _day_flags.filter(pl.col("has_overnight_data")).height
    _n_imv = _day_flags.filter(pl.col("has_overnight_data") & pl.col("has_imv")).height

    consort_partial = {
        "total": _total, "n_overnight": _n_overnight, "n_imv": _n_imv,
    }
    return consort_partial, failure_flags_df, failure_reason_df


@app.cell
def _(np, pl, tqdm):
    def process_cohort(df: pl.DataFrame) -> pl.DataFrame:
        """Find patient-days with >= 6h continuous SBT eligibility (hemodynamic + respiratory stability) in the 22:00-06:00 window."""
        df = df.sort("hospitalization_id", "recorded_dttm")

        # Mark each row: IMV + no paralytics + hemodynamic + respiratory stability
        df = df.with_columns(
            (
                (pl.col("device_category").str.to_lowercase() == "imv")
                & (pl.col("max_paralytics") <= 0)
                & (pl.col("hemodynamic_stable") == 1)
                & (pl.col("respiratory_stable") == 1)
            ).cast(pl.Int32).alias("all_conditions_check")
        )

        # Only keep days that have at least one IMV observation
        vented_keys = (
            df.filter(pl.col("device_category").str.to_lowercase() == "imv")
            .select("hosp_id_icu_day")
            .unique()
        )
        df = df.join(vented_keys, on="hosp_id_icu_day", how="semi")

        _6h_us = int(6 * 3600 * 1e6)  # 6 hours in microseconds
        _22h = np.timedelta64(22, "h")
        _1d = np.timedelta64(1, "D")
        _6h = np.timedelta64(6, "h")

        result_hosp = []
        result_day = []
        result_time = []

        # Loop over each patient
        for (hosp_id,), hosp_df in tqdm(df.group_by("hospitalization_id"), desc="SBT eligibility check (hemo+resp)"):
            hosp_df = hosp_df.sort("recorded_dttm")
            times = hosp_df["recorded_dttm"].to_numpy().astype("datetime64[us]")
            conditions = hosp_df["all_conditions_check"].to_numpy()
            vent_days = hosp_df["icu_day_date"].unique().sort().to_numpy().astype("datetime64[us]")

            # For each vent-day, look at the overnight window (22:00 → 06:00)
            for date in vent_days:
                # Build the 8-hour overnight window
                start_time = date + _22h - _1d   # previous day 22:00
                end_time = date + _6h             # current day 06:00
                mask = (times >= start_time) & (times <= end_time)
                w_times = times[mask]
                w_conds = conditions[mask]

                # Skip if no data in window or none of the conditions are ever met
                if len(w_times) == 0 or not w_conds.any():
                    continue

                # Identify contiguous segments where conditions flip on/off
                changes = np.concatenate([[True], np.diff(w_conds) != 0])
                groups = np.cumsum(changes)

                # Look only at segments where all conditions ARE met (==1)
                for g in np.unique(groups[w_conds == 1]):
                    seg_mask = groups == g
                    seg_times = w_times[seg_mask]

                    if len(seg_times) < 2:
                        continue

                    # Cumulative duration from first observation in segment
                    seg_i = seg_times.astype("datetime64[us]").astype(np.int64)
                    cum_dur = seg_i - seg_i[0]

                    # If this segment lasts >= 6 hours, the day is eligible
                    if cum_dur[-1] >= _6h_us:
                        # Record the exact timestamp when we hit 6 hours
                        hit_idx = np.searchsorted(cum_dur, _6h_us)
                        result_hosp.append(hosp_id)
                        result_day.append(date)
                        result_time.append(seg_times[hit_idx])
                        break  # one qualifying segment is enough for this day

        if not result_hosp:
            return pl.DataFrame(schema={
                "hospitalization_id": df.schema["hospitalization_id"],
                "current_day_key": df.schema["icu_day_date"],
                "event_time_at_6_hours": df.schema["recorded_dttm"],
            })
        return pl.DataFrame({
            "hospitalization_id": result_hosp,
            "current_day_key": np.array(result_day, dtype="datetime64[us]"),
            "event_time_at_6_hours": np.array(result_time, dtype="datetime64[us]"),
        })

    return (process_cohort,)


@app.cell
def _(df, failure_reason_df, pl, process_cohort):
    # Run the eligibility check
    result_df = process_cohort(df)
    print(f"Encounter-days with ≥6h eligibility: {result_df.height:,}")

    # Left-join eligibility info onto every row (non-eligible days get NaN)
    cohort_elig = df.join(
        result_df.select("hospitalization_id", "current_day_key", "event_time_at_6_hours"),
        left_on=["hospitalization_id", "icu_day_date"],
        right_on=["hospitalization_id", "current_day_key"],
        how="left",
    )

    # Find the FIRST observation on each eligible day that occurs AFTER the 6h mark
    _first_times = (
        cohort_elig
        .filter(
            pl.col("event_time_at_6_hours").is_not_null()
            & (pl.col("recorded_dttm") >= pl.col("event_time_at_6_hours"))
        )
        .group_by("hospitalization_id", "icu_day_date")
        .agg(pl.col("recorded_dttm").min().alias("_first_eligible_time"))
    )

    # Flag just that first post-threshold row
    cohort_elig = cohort_elig.join(_first_times, on=["hospitalization_id", "icu_day_date"], how="left")
    cohort_elig = cohort_elig.with_columns(
        pl.when(pl.col("recorded_dttm") == pl.col("_first_eligible_time"))
        .then(1.0)
        .otherwise(None)
        .alias("eligible_event")
    )

    # Remove flag from the very last observation per patient
    # (can't assess an SBT if there's no follow-up data after it)
    cohort_elig = cohort_elig.with_columns(
        pl.when(
            pl.col("recorded_dttm") == pl.col("recorded_dttm").max().over("hospitalization_id")
        )
        .then(None)
        .otherwise(pl.col("eligible_event"))
        .alias("eligible_event")
    )

    # Mark ALL rows on eligible days so we can filter to them later
    _eligible_keys = (
        cohort_elig.filter(pl.col("eligible_event") == 1)
        .select("hosp_id_icu_day")
        .unique()
    )
    cohort_elig = cohort_elig.with_columns(
        pl.col("hosp_id_icu_day").is_in(_eligible_keys["hosp_id_icu_day"].to_list())
        .cast(pl.Int32)
        .alias("on_vent_eligible")
    )
    cohort_elig = cohort_elig.drop("event_time_at_6_hours", "_first_eligible_time")

    # Join eligibility failure reason; null out for eligible days
    cohort_elig = cohort_elig.join(failure_reason_df, on="hosp_id_icu_day", how="left")
    cohort_elig = cohort_elig.with_columns(
        pl.when(pl.col("on_vent_eligible") == 1)
        .then(pl.lit(None).cast(pl.Utf8))
        .otherwise(pl.col("eligibility_failure_reason"))
        .alias("eligibility_failure_reason")
    )

    n_eligible_days = cohort_elig.filter(pl.col("on_vent_eligible") == 1)["hosp_id_icu_day"].n_unique()
    print(f"Eligible hosp-day keys: {n_eligible_days:,}")
    return cohort_elig, n_eligible_days


@app.cell
def _(cohort_elig, np, pl, tqdm):
    # ── SBT Delivery: Ventilator Mode FLIP Detection ──
    # A FLIP = ventilator mode switches from controlled to PS/CPAP/T-piece
    # with pressure_support_set <= 8 and peep_set <= 8
    _FLAG_COLS = ["SBT_delivery_2min_primary", "SBT_delivery_5min_secondary"]

    # Keep only rows from eligible days
    vent_eligible = (
        cohort_elig.filter(pl.col("on_vent_eligible") == 1)
        .sort("hospitalization_id", "recorded_dttm")
    )

    # Re-join the eligibility threshold time for FLIP timing check
    # (FLIP must occur AFTER the 6h threshold was reached)
    _elig_times = (
        cohort_elig.filter(pl.col("eligible_event") == 1)
        .select("hosp_id_icu_day", "recorded_dttm")
        .rename({"recorded_dttm": "_elig_threshold_time"})
    )

    vent_eligible = vent_eligible.join(_elig_times, on="hosp_id_icu_day", how="left")

    # ── Precompute FLIP condition vectorized ──
    # FLIP = device is IMV AND (
    #   (mode_category contains "pressure support" or "cpap"
    #    AND pressure_support_set <= 8 AND peep_set <= 8)
    #   OR mode_name matches t-piece pattern
    # )
    vent_eligible = vent_eligible.with_columns(
        (
            (pl.col("device_category").str.to_lowercase() == "imv")
            & (
                (
                    (
                        pl.col("mode_category").fill_null("").str.to_lowercase().str.contains("pressure support")
                        | pl.col("mode_category").fill_null("").str.to_lowercase().str.contains("cpap")
                    )
                    & (pl.col("pressure_support_set").fill_null(float("inf")) <= 8)
                    & (pl.col("peep_set").fill_null(float("inf")) <= 8)
                )
                | pl.col("mode_name").fill_null("").str.to_lowercase().str.contains(r"^t-?piece$")
            )
        ).alias("flip_check_flag")
    )

    # Add row index for joining flag results back
    vent_eligible = vent_eligible.with_row_index("_row_idx")

    # ── Main loop: evaluate FLIP at each candidate timepoint ──
    _delta_2min_us = int(2 * 60 * 1e6)    # 2 min in microseconds
    _delta_5min_us = int(5 * 60 * 1e6)    # 5 min in microseconds
    _flag_results = []

    for _key in tqdm(vent_eligible["hosp_id_icu_day"].unique().sort().to_list(), desc="Evaluating SBT FLIP"):
        _grp = vent_eligible.filter(pl.col("hosp_id_icu_day") == _key).sort("recorded_dttm")
        _n = _grp.height
        if _n == 0:
            continue

        _row_idxs = _grp["_row_idx"].to_numpy()
        _times = _grp["recorded_dttm"].to_numpy().astype("datetime64[us]").astype(np.int64)
        _flips = _grp["flip_check_flag"].to_numpy()

        # Get the eligibility threshold time for this day
        _thresh_vals = _grp["_elig_threshold_time"].to_numpy().astype("datetime64[us]").astype(np.int64)
        _thresh_time = _thresh_vals[0]  # same for all rows in this day

        _found_2min = False
        _found_5min = False

        # Iterate through candidate FLIP timepoints (must be after 6h threshold)
        for i in range(_n):
            if _found_2min and _found_5min:
                break

            if not _flips[i]:
                continue

            # FLIP must occur AFTER the eligibility threshold time
            if _times[i] <= _thresh_time:
                continue

            cur_time = _times[i]

            # Check 2-min sustained FLIP
            if not _found_2min:
                fw_2min_mask = (_times >= cur_time) & (_times <= cur_time + _delta_2min_us)
                fw_flips = _flips[fw_2min_mask]
                if len(fw_flips) > 0 and fw_flips.all():
                    _found_2min = True
                    _row = {"_row_idx": int(_row_idxs[i]), "SBT_delivery_2min_primary": 1.0}
                    # Check if this same point also satisfies 5min
                    if not _found_5min:
                        fw_5min_mask = (_times >= cur_time) & (_times <= cur_time + _delta_5min_us)
                        fw_flips_5 = _flips[fw_5min_mask]
                        if len(fw_flips_5) > 0 and fw_flips_5.all():
                            _found_5min = True
                            _row["SBT_delivery_5min_secondary"] = 1.0
                    _flag_results.append(_row)
                    continue

            # Check 5-min sustained FLIP (independent search)
            if not _found_5min:
                fw_5min_mask = (_times >= cur_time) & (_times <= cur_time + _delta_5min_us)
                fw_flips_5 = _flips[fw_5min_mask]
                if len(fw_flips_5) > 0 and fw_flips_5.all():
                    _found_5min = True
                    _flag_results.append({
                        "_row_idx": int(_row_idxs[i]),
                        "SBT_delivery_5min_secondary": 1.0,
                    })

    # Build flag result DataFrame and join back
    if _flag_results:
        _flag_df = pl.DataFrame(
            _flag_results,
            schema={"_row_idx": pl.UInt32, **{c: pl.Float64 for c in _FLAG_COLS}},
        )
        # Multiple rows may exist for the same day (2min and 30min at different times)
        # Group by _row_idx and take max to avoid duplicates
        _flag_df = _flag_df.group_by("_row_idx").agg(
            [pl.col(c).max().alias(c) for c in _FLAG_COLS]
        )
        vent_eligible = vent_eligible.join(_flag_df, on="_row_idx", how="left")
    else:
        vent_eligible = vent_eligible.with_columns(
            *[pl.lit(None).cast(pl.Float64).alias(_f) for _f in _FLAG_COLS]
        )

    vent_eligible = vent_eligible.drop("_row_idx", "_elig_threshold_time")

    sbt_flag_cols = _FLAG_COLS
    return sbt_flag_cols, vent_eligible


@app.cell
def _(DATA_DIR, Path, TIMEZONE, cohort, cohort_elig, df, failure_flags_df, pl, sbt_flag_cols, vent_eligible):
    cohort_with_delivery = vent_eligible

    # Aggregate: for each hosp_id_icu_day, take the MAX of each flag column
    # (so if any single row in a day had flag=1, the whole day gets flag=1)

    # Step 1: ALL days — base columns from cohort_elig
    _gt_cols = ["sbt_screen_pass_fail", "sbt_screen_performed", "sbt_delivery_pass_fail", "sbt_delivery_performed"]
    _base_cols = _gt_cols + ["eligible_event"]
    _all_days = (
        cohort_elig.select(["hosp_id_icu_day"] + _base_cols)
        .group_by("hosp_id_icu_day")
        .agg([pl.col(c).max().alias(c) for c in _base_cols])
    )

    # --- Vent-day filter ---
    # Condition 1: day is on or after intubation (calendar date)
    _day_info = (
        cohort_elig.select("hosp_id_icu_day", "hospitalization_id", "icu_day_date")
        .unique(subset=["hosp_id_icu_day"])
        .join(
            cohort.select("hospitalization_id", "intubation_time"),
            on="hospitalization_id",
            how="left",
        )
    )
    _post_intub_keys = (
        _day_info.filter(
            pl.col("intubation_time").dt.truncate("1d") <= pl.col("icu_day_date")
        )
        .select("hosp_id_icu_day")
    )

    # Condition 2: IMV charted in overnight window (prev day 22:00 – today 06:00)
    _imv_obs = (
        cohort_elig.filter(pl.col("device_category").str.to_lowercase() == "imv")
        .select("hospitalization_id", "recorded_dttm")
    )
    _day_keys = (
        cohort_elig.select("hosp_id_icu_day", "hospitalization_id", "icu_day_date")
        .unique(subset=["hosp_id_icu_day"])
    )
    _imv_day_keys = (
        _day_keys.join(_imv_obs, on="hospitalization_id", how="inner")
        .filter(
            (pl.col("recorded_dttm") >= (pl.col("icu_day_date") - pl.duration(hours=2)))
            & (pl.col("recorded_dttm") <= (pl.col("icu_day_date") + pl.duration(hours=6)))
        )
        .select("hosp_id_icu_day").unique()
    )

    # Vent day = both conditions met; keep ALL ICU days with a flag
    _vent_day_keys = _post_intub_keys.join(_imv_day_keys, on="hosp_id_icu_day", how="semi")
    _all_days = _all_days.join(
        _vent_day_keys.with_columns(pl.lit(True).alias("is_vent_day")),
        on="hosp_id_icu_day", how="left",
    ).with_columns(pl.col("is_vent_day").fill_null(False))

    # Step 2: Eligible days only — SBT flag columns from cohort_with_delivery
    _flag_agg = (
        cohort_with_delivery.select(["hosp_id_icu_day"] + sbt_flag_cols)
        .group_by("hosp_id_icu_day")
        .agg([pl.col(c).max().alias(c) for c in sbt_flag_cols])
    )

    # Step 3: Left-join flags onto all days
    df_grouped = _all_days.join(_flag_agg, on="hosp_id_icu_day", how="left").sort("hosp_id_icu_day")

    df_grouped = df_grouped.with_columns(
        pl.when(pl.col("eligible_event") == 1)
        .then(True)
        .otherwise(False)
        .alias("is_eligible")
    )

    # Join independent failure flags + mark no_6h_window
    df_grouped = df_grouped.join(failure_flags_df, on="hosp_id_icu_day", how="left")
    df_grouped = df_grouped.with_columns(
        (~pl.col("is_eligible")).alias("no_6h_window"),
    )

    # Build the ground truth: a day counts as "delivered" if ANY of the 4 EHR
    # flowsheet columns is charted (non-null), AND the day was eligible
    df_grouped = df_grouped.with_columns(
        pl.when(
            (
                pl.col("sbt_screen_pass_fail").is_not_null()
                | pl.col("sbt_screen_performed").is_not_null()
                | pl.col("sbt_delivery_pass_fail").is_not_null()
                | pl.col("sbt_delivery_performed").is_not_null()
            )
            & (pl.col("eligible_event") == 1)
        )
        .then(1.0)
        .otherwise(None)
        .alias("sbt_ground_truth")
    )

    # ── Delivery failure reasons ──
    # Did a FLIP ever occur on this day (after threshold)?
    _flip_occurred = (
        vent_eligible.group_by("hosp_id_icu_day").agg(
            pl.col("flip_check_flag").any().alias("_has_flip"),
        )
    )

    df_grouped = df_grouped.join(_flip_occurred, on="hosp_id_icu_day", how="left")
    df_grouped = df_grouped.with_columns(
        # 2-min delivery failure
        pl.when(pl.col("SBT_delivery_2min_primary") == 1).then(None)
        .when(pl.col("eligible_event").is_null()).then(None)
        .when(pl.col("_has_flip").not_()).then(pl.lit("No mode FLIP to PS/CPAP/T-piece"))
        .otherwise(pl.lit("FLIP not sustained for 2min"))
        .alias("SBT_delivery_2min_primary_failure"),
        # 5-min delivery failure
        pl.when(pl.col("SBT_delivery_5min_secondary") == 1).then(None)
        .when(pl.col("eligible_event").is_null()).then(None)
        .when(pl.col("_has_flip").not_()).then(pl.lit("No mode FLIP to PS/CPAP/T-piece"))
        .otherwise(pl.lit("FLIP not sustained for 5min"))
        .alias("SBT_delivery_5min_secondary_failure"),
    ).drop("_has_flip")

    # ── Join eligibility failure reason (one per day) ──
    _elig_reason = (
        cohort_elig.select("hosp_id_icu_day", "eligibility_failure_reason")
        .group_by("hosp_id_icu_day").agg(pl.col("eligibility_failure_reason").first())
    )
    df_grouped = df_grouped.join(_elig_reason, on="hosp_id_icu_day", how="left")

    # ── Keep all days; add eligibility pass flag ──
    df_grouped = df_grouped.with_columns(
        pl.col("is_eligible").alias("eligibility_day_pass"),
    )

    # Set failure reason for non-vent days
    df_grouped = df_grouped.with_columns(
        pl.when(~pl.col("is_vent_day") & pl.col("eligibility_failure_reason").is_null())
        .then(pl.lit("Not a vent day"))
        .otherwise(pl.col("eligibility_failure_reason"))
        .alias("eligibility_failure_reason")
    )

    # Capture total vent day count (vent days only)
    n_vent_days = df_grouped.filter(pl.col("is_vent_day")).height

    # ── Enrich with demographics, priors, and outcome ──
    from clifpy.utils.comorbidity import calculate_cci
    from clifpy.utils.sofa_polars import compute_sofa_polars

    # Map hosp_id_icu_day → hospitalization_id + icu_day + icu_day_date
    _key_map = (
        df.select("hosp_id_icu_day", "hospitalization_id", "icu_day", "icu_day_date")
        .unique(subset=["hosp_id_icu_day"])
    )
    df_grouped = df_grouped.join(_key_map, on="hosp_id_icu_day", how="left")

    # Sequential vent_day numbering (1, 2, 3... for vent days only; null for non-vent)
    _vent_day_ranking = (
        df_grouped.filter(pl.col("is_vent_day"))
        .select("hosp_id_icu_day", "hospitalization_id", "icu_day")
        .sort("hospitalization_id", "icu_day")
        .with_columns(
            pl.col("icu_day")
            .rank("dense")
            .over("hospitalization_id")
            .cast(pl.Int64)
            .alias("vent_day")
        )
        .select("hosp_id_icu_day", "vent_day")
    )
    df_grouped = df_grouped.join(_vent_day_ranking, on="hosp_id_icu_day", how="left")

    # Charlson Comorbidity Index (ICD-10 only)
    _dx = pl.read_parquet(Path(DATA_DIR) / "clif_hospital_diagnosis.parquet")
    _dx = _dx.filter(
        pl.col("diagnosis_code_format") == "ICD10CM",
        pl.col("hospitalization_id").is_in(cohort["hospitalization_id"].to_list()),
    )
    _cci = pl.from_pandas(
        calculate_cci(_dx, hierarchy=True)
    ).select("hospitalization_id", "cci_score")

    # SOFA — worst per vent-day, then shift to prior day
    _day_cohort = (
        df.select("hosp_id_icu_day", "hospitalization_id", "icu_day_date")
        .unique(subset=["hosp_id_icu_day"])
        .with_columns(
            pl.col("icu_day_date").alias("start_dttm"),
            (pl.col("icu_day_date") + pl.duration(days=1)).alias("end_dttm"),
        )
    )
    _daily_sofa = compute_sofa_polars(
        data_directory=DATA_DIR,
        cohort_df=_day_cohort,
        filetype="parquet",
        id_name="hosp_id_icu_day",
        extremal_type="worst",
        fill_na_scores_with_zero=True,
        remove_outliers=True,
        timezone=TIMEZONE,
    ).select("hosp_id_icu_day", "sofa_total")

    _sofa_with_day = _daily_sofa.join(
        df.select("hosp_id_icu_day", "hospitalization_id", "icu_day")
        .unique(subset=["hosp_id_icu_day"]),
        on="hosp_id_icu_day", how="left",
    )
    _sofa_prior = _sofa_with_day.with_columns(
        (pl.col("icu_day") + 1).alias("_next")
    ).select("hospitalization_id", "_next", pl.col("sofa_total").alias("sofa_prior"))

    # Prior-day aggregates from wide dataset
    _day_agg = df.group_by(
        "hosp_id_icu_day", "hospitalization_id", "icu_day"
    ).agg(
        ((pl.col("lorazepam") > 0) | (pl.col("midazolam") > 0))
        .any()
        .alias("bzd_prior"),
        (pl.col("propofol") > 0).any().alias("propofol_prior"),
        (
            (pl.col("fentanyl") > 0)
            | (pl.col("hydromorphone") > 0)
            | (pl.col("morphine") > 0)
        )
        .any()
        .alias("opioid_prior"),
        (pl.col("max_paralytics") > 0).any().alias("nmb_prior"),
        # Individual sedation med flags
        (pl.col("fentanyl") > 0).any().alias("fentanyl_prior"),
        (pl.col("hydromorphone") > 0).any().alias("hydromorphone_prior"),
        (pl.col("morphine") > 0).any().alias("morphine_prior"),
        (pl.col("lorazepam") > 0).any().alias("lorazepam_prior"),
        (pl.col("midazolam") > 0).any().alias("midazolam_prior"),
        pl.col("rass").median().alias("rass_prior"),
        pl.col("nee").median().alias("nee_prior"),
        pl.col("fio2_set").median().alias("fio2_prior"),
        pl.col("peep_set").median().alias("peep_prior"),
        pl.col("gcs_total").median().alias("gcs_prior"),
    )
    _prior_cols = [
        "bzd_prior", "propofol_prior", "opioid_prior", "nmb_prior",
        "fentanyl_prior", "hydromorphone_prior", "morphine_prior",
        "lorazepam_prior", "midazolam_prior",
        "rass_prior", "nee_prior", "fio2_prior", "peep_prior", "gcs_prior",
    ]
    _prior = _day_agg.with_columns(
        (pl.col("icu_day") + 1).alias("_next")
    ).select("hospitalization_id", "_next", *_prior_cols)
    df_grouped = df_grouped.join(
        _prior,
        left_on=["hospitalization_id", "icu_day"],
        right_on=["hospitalization_id", "_next"],
        how="left",
    )

    # Demographics, CCI, SOFA
    _demo = cohort.select(
        "hospitalization_id", "age_at_admission", "sex_category",
        "race_category", "ethnicity_category", "language_category",
        "weight_kg", "height_cm", "bmi",
        "hospital_los_hours", "first_icu_los_hours", "imv_duration_hours",
        "discharge_category",
    )
    df_grouped = df_grouped.join(_demo, on="hospitalization_id", how="left")
    df_grouped = df_grouped.join(_cci, on="hospitalization_id", how="left")
    df_grouped = df_grouped.join(
        _sofa_prior,
        left_on=["hospitalization_id", "icu_day"],
        right_on=["hospitalization_id", "_next"],
        how="left",
    )
    df_grouped = df_grouped.with_columns(
        (pl.col("sex_category").str.to_lowercase() == "female").alias("is_female"),
        (pl.col("language_category").str.to_lowercase() == "english").alias("is_english"),
    )

    # Outcome: 3=death, 2=extubated, 0=alive (per day)
    _patient = pl.read_parquet(Path(DATA_DIR) / "clif_patient.parquet").select(
        "patient_id", "death_dttm"
    )
    _death_dtype = _patient.schema["death_dttm"]
    if _death_dtype == pl.Date:
        _patient = _patient.with_columns(
            pl.col("death_dttm").cast(pl.Datetime("us"))
        )
    elif isinstance(_death_dtype, pl.Datetime) and _death_dtype.time_zone is not None:
        _patient = _patient.with_columns(
            pl.col("death_dttm").dt.replace_time_zone(None)
        )
    _death_info = cohort.select(
        "hospitalization_id", "patient_id", "extubation_time", "discharge_dttm",
    ).join(_patient, on="patient_id", how="left")
    df_grouped = df_grouped.join(
        _death_info.select("hospitalization_id", "extubation_time", "discharge_dttm", "death_dttm"),
        on="hospitalization_id", how="left",
    )
    df_grouped = df_grouped.with_columns(
        pl.when(
            pl.col("death_dttm").is_not_null()
            & (pl.col("death_dttm").dt.truncate("1d") == pl.col("icu_day_date"))
        ).then(3)
        .when(
            (pl.col("discharge_category").str.to_lowercase() == "expired")
            & pl.col("death_dttm").is_null()
            & (pl.col("discharge_dttm").dt.truncate("1d") == pl.col("icu_day_date"))
        ).then(3)
        .when(
            pl.col("extubation_time").is_not_null()
            & (pl.col("extubation_time").dt.truncate("1d") == pl.col("icu_day_date"))
        ).then(2)
        .otherwise(0)
        .alias("outcome")
    )

    print(f"Day-level aggregated: {df_grouped.height:,} hosp-days")
    for _col in sbt_flag_cols + ["sbt_ground_truth"]:
        _n = df_grouped.filter(pl.col(_col) == 1).height
        print(f"  {_col}: {_n:,} days flagged")
    return cohort_with_delivery, df_grouped, n_vent_days


@app.cell
def _(
    OUTPUT_SHARE,
    SITE_NAME,
    consort_partial,
    df_grouped,
    json,
    pl,
    plot_consort,
    plot_upset,
):
    _n_vent_days = df_grouped.filter(pl.col("is_vent_day")).height
    _n_elig = df_grouped.filter(pl.col("eligible_event") == 1).height
    _n_not_elig = _n_vent_days - _n_elig

    consort_sbt = [
        {"step": 0, "description": "Total patient days in index ICU",
         "remaining": consort_partial["total"], "excluded": 0, "reason": None},
        {"step": 1, "description": "Patient-days with data (22:00–06:00)",
         "remaining": consort_partial["n_overnight"],
         "excluded": consort_partial["total"] - consort_partial["n_overnight"],
         "reason": "No overnight data"},
        {"step": 2, "description": "Vent-days with IMV",
         "remaining": _n_vent_days,
         "excluded": consort_partial["n_overnight"] - _n_vent_days,
         "reason": "No IMV"},
        {"step": 3, "description": "Eligible days (≥6h IMV, no paralytics, both stable)",
         "remaining": _n_elig,
         "excluded": _n_not_elig, "reason": "Not eligible"},
    ]

    with open(OUTPUT_SHARE / f"consort_sbt_{SITE_NAME}.json", "w") as _f:
        json.dump(consort_sbt, _f, indent=2)
    plot_consort(consort_sbt, OUTPUT_SHARE / f"consort_sbt_{SITE_NAME}.png")

    # Delivery sub-analysis print
    _n_2min = df_grouped.filter(pl.col("SBT_delivery_2min_primary") == 1).height
    _n_5min = df_grouped.filter(pl.col("SBT_delivery_5min_secondary") == 1).height
    print(f"SBT CONSORT saved. Eligible days: {_n_elig:,}")
    print(f"  2min delivery: {_n_2min:,}  |  5min delivery: {_n_5min:,}")

    # Eligibility failure CSV + UpSet plot
    _fail_cols = ["has_paralytics", "resp_unstable", "hemo_unstable", "no_6h_window"]
    _fail_df = df_grouped.filter(
        ~pl.col("is_eligible")
        & ~pl.col("eligibility_failure_reason").is_in(["No overnight data", "No IMV"])
    ).select(["hosp_id_icu_day"] + _fail_cols)
    _fail_df.write_csv(OUTPUT_SBT / f"eligibility_failures_{SITE_NAME}.csv")
    plot_upset(_fail_df.to_pandas(), _fail_cols, OUTPUT_SHARE / f"upset_eligibility_{SITE_NAME}.png")
    print(f"Eligibility failure CSV + UpSet plot saved ({_fail_df.height:,} non-eligible days)")

    # Detection failure UpSet (eligible days where SBT not detected)
    _det_df = (
        df_grouped.filter(pl.col("is_eligible"))
        .filter(pl.col("SBT_delivery_5min_secondary").is_null())
        .with_columns(
            (pl.col("SBT_delivery_2min_primary_failure") == "No mode FLIP to PS/CPAP/T-piece").fill_null(False).alias("no_flip"),
            (pl.col("SBT_delivery_2min_primary_failure") == "FLIP not sustained for 2min").fill_null(False).alias("flip_not_sustained_2min"),
            (pl.col("SBT_delivery_2min_primary").is_not_null() & pl.col("SBT_delivery_5min_secondary").is_null()).alias("flip_sustained_2min_not_5min"),
        )
    )
    _det_cols = ["no_flip", "flip_not_sustained_2min", "flip_sustained_2min_not_5min"]
    _det_df.select(["hosp_id_icu_day"] + _det_cols).write_csv(OUTPUT_SBT / f"detection_failures_{SITE_NAME}.csv")
    plot_upset(_det_df.select(_det_cols).to_pandas(), _det_cols, OUTPUT_SHARE / f"upset_detection_{SITE_NAME}.png")
    print(f"Detection failure UpSet saved ({_det_df.height:,} eligible days with failed detection)")
    return


@app.cell
def _(OUTPUT_SHARE, df_grouped, np, pd, pl, sbt_flag_cols):
    from sklearn.metrics import cohen_kappa_score, confusion_matrix
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams["font.family"] = "Arial"
    rcParams["font.size"] = 8

    def _bootstrap_metrics_ci(y_true, y_pred, n_boot=2000, ci=0.95, seed=42):
        """Bootstrap 95% CIs for all concordance metrics."""
        rng = np.random.default_rng(seed)
        n = len(y_true)
        results = {k: [] for k in [
            "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "F1", "Cohen_Kappa"
        ]}
        for _ in range(n_boot):
            idx = rng.choice(n, size=n, replace=True)
            yt, yp = y_true.iloc[idx], y_pred.iloc[idx]
            cm = confusion_matrix(yt, yp, labels=[0, 1])
            tn, fp, fn, tp = cm.ravel()
            tot = cm.sum()
            results["Accuracy"].append((tp + tn) / tot)
            results["Sensitivity"].append(tp / (tp + fn) if tp + fn else np.nan)
            results["Specificity"].append(tn / (tn + fp) if tn + fp else np.nan)
            results["PPV"].append(tp / (tp + fp) if tp + fp else np.nan)
            results["NPV"].append(tn / (tn + fn) if tn + fn else np.nan)
            p = tp / (tp + fp) if tp + fp else 0
            r = tp / (tp + fn) if tp + fn else 0
            results["F1"].append(2 * p * r / (p + r) if p + r else np.nan)
            try:
                results["Cohen_Kappa"].append(cohen_kappa_score(yt, yp))
            except Exception:
                results["Cohen_Kappa"].append(np.nan)
        alpha = (1 - ci) / 2
        ci_dict = {}
        for k, vals in results.items():
            arr = np.array([v for v in vals if not np.isnan(v)])
            if len(arr):
                ci_dict[k] = (float(np.percentile(arr, alpha * 100)),
                              float(np.percentile(arr, (1 - alpha) * 100)))
            else:
                ci_dict[k] = (np.nan, np.nan)
        return ci_dict

    def _landis_koch(kappa):
        """Classify kappa value using the Landis-Koch agreement scale."""
        if kappa < 0: return "Poor"
        elif kappa < 0.21: return "Slight"
        elif kappa < 0.41: return "Fair"
        elif kappa < 0.61: return "Moderate"
        elif kappa < 0.81: return "Substantial"
        return "Almost Perfect"

    # Fill null flags with 0 in Polars, then convert to pandas once for sklearn/matplotlib
    _fill_cols = sbt_flag_cols + ["sbt_ground_truth"]
    _con_df = (
        df_grouped
        .filter(pl.col("is_vent_day") & (pl.col("eligible_event") == 1))
        .with_columns([pl.col(c).fill_null(0) for c in _fill_cols])
        .to_pandas()
    )

    _metrics_wide = {}

    # Compare each algorithmic flag against the flowsheet ground truth
    for _col in sbt_flag_cols:
        _y_true = _con_df["sbt_ground_truth"]  # ground truth (from EHR flowsheet)
        _y_pred = _con_df[_col]                 # algorithmic prediction

        # Compute confusion matrix and standard metrics
        _cm = confusion_matrix(_y_true, _y_pred, labels=[0, 1])
        _tn, _fp, _fn, _tp = _cm.ravel()
        _total = _cm.sum()
        _accuracy = (_tp + _tn) / _total
        _precision = _tp / (_tp + _fp) if _tp + _fp else 0
        _recall = _tp / (_tp + _fn) if _tp + _fn else 0
        _f1 = 2 * _precision * _recall / (_precision + _recall) if _precision + _recall else 0
        _specificity = _tn / (_tn + _fp) if _tn + _fp else 0
        _npv = _tn / (_tn + _fn) if _tn + _fn else 0
        _kappa = cohen_kappa_score(_y_true, _y_pred)

        # Bootstrap 95% CIs for all metrics
        _ci = _bootstrap_metrics_ci(_y_true, _y_pred)

        # Save labelled confusion matrix CSV
        _cm_df = pd.DataFrame(
            _cm,
            index=["Actual Negative", "Actual Positive"],
            columns=["Predicted Negative", "Predicted Positive"],
        )
        _cm_df.to_csv(OUTPUT_SHARE / f"confusion_matrix_{_col}.csv")

        # Save a JAMA-style confusion matrix heatmap
        _cm_pct = _cm / _total * 100
        _fig, _ax = plt.subplots(figsize=(3.5, 3.0))
        _ax.imshow(_cm, cmap="cividis")
        _ax.set_xticks([0, 1]); _ax.set_yticks([0, 1])
        _ax.set_xticklabels(["No Delivery", "Delivery"])
        _ax.set_yticklabels(["No Delivery", "Delivery"])
        _ax.set_xlabel(f"{_col} flag"); _ax.set_ylabel("Flowsheet delivery flag")
        _ax.set_title(f"Concordance: flowsheet vs {_col}", fontweight="bold")
        _ax.spines["top"].set_visible(False); _ax.spines["right"].set_visible(False)
        for _i in range(2):
            for _j in range(2):
                _ax.text(_j, _i, f"{_cm[_i, _j]}\n({_cm_pct[_i, _j]:.1f}%)",
                         ha="center", va="center",
                         color="white" if _cm[_i, _j] > _cm.max() / 2 else "black")
        _fig.tight_layout()
        _fig.savefig(OUTPUT_SHARE / f"confusion_matrix_{_col}.png", bbox_inches="tight", dpi=300)
        plt.close(_fig)

        # Build row-wide metrics column for this flag
        _metrics_wide[_col] = {
            "TP": _tp,
            "FP": _fp,
            "FN": _fn,
            "TN": _tn,
            "Accuracy": f"{_accuracy:.4f} ({_ci['Accuracy'][0]:.4f}\u2013{_ci['Accuracy'][1]:.4f})",
            "Sensitivity (Recall)": f"{_recall:.4f} ({_ci['Sensitivity'][0]:.4f}\u2013{_ci['Sensitivity'][1]:.4f})",
            "Specificity": f"{_specificity:.4f} ({_ci['Specificity'][0]:.4f}\u2013{_ci['Specificity'][1]:.4f})",
            "PPV (Precision)": f"{_precision:.4f} ({_ci['PPV'][0]:.4f}\u2013{_ci['PPV'][1]:.4f})",
            "NPV": f"{_npv:.4f} ({_ci['NPV'][0]:.4f}\u2013{_ci['NPV'][1]:.4f})",
            "F1": f"{_f1:.4f} ({_ci['F1'][0]:.4f}\u2013{_ci['F1'][1]:.4f})",
            "Cohen_Kappa": f"{_kappa:.4f} ({_ci['Cohen_Kappa'][0]:.4f}\u2013{_ci['Cohen_Kappa'][1]:.4f})",
            "Kappa_Interpretation": _landis_koch(_kappa),
        }

    # Save row-wide summary (metrics as rows, flags as columns)
    concordance_df = pd.DataFrame(_metrics_wide)
    concordance_df.index.name = "Metric"
    concordance_df.to_csv(OUTPUT_SHARE / "delivery_concordance_summary.csv")
    print(concordance_df.to_string())
    return


@app.cell
def _(
    OUTPUT_SHARE,
    SITE_NAME,
    df_grouped,
    json,
    pd,
    pl,
):
    _tbl = df_grouped

    # ── Formatting helpers ──
    def _median_iqr(s):
        s = s.drop_nulls().cast(pl.Float64)
        if len(s) == 0:
            return "—"
        return f"{s.median():.1f} ({s.quantile(0.25):.1f}–{s.quantile(0.75):.1f})"

    def _mean_sd(s):
        s = s.drop_nulls().cast(pl.Float64)
        if len(s) == 0:
            return "—"
        return f"{s.mean():.1f} ({s.std():.1f})"

    def _n_pct(mask, n_total):
        n = int(mask.drop_nulls().sum())
        return f"{n} ({n / n_total * 100:.1f}%)" if n_total else "—"

    def _cat_rows(s, n_total):
        s = s.drop_nulls()
        if len(s) == 0:
            return []
        vc = s.value_counts().sort("count", descending=True)
        col = [c for c in vc.columns if c != "count"][0]
        return [
            (f"  {r[col]}", f"{r['count']} ({r['count'] / n_total * 100:.1f}%)")
            for r in vc.iter_rows(named=True)
        ]

    # ── Step 6b: Raw-value helpers for machine-parseable JSON ──
    def _median_iqr_raw(s):
        n_missing = int(s.is_null().sum())
        s = s.drop_nulls().cast(pl.Float64)
        if len(s) == 0:
            return {"median": None, "q25": None, "q75": None, "n_missing": n_missing}
        return {
            "median": round(float(s.median()), 1),
            "q25": round(float(s.quantile(0.25)), 1),
            "q75": round(float(s.quantile(0.75)), 1),
            "n_missing": n_missing,
        }

    def _mean_sd_raw(s):
        n_missing = int(s.is_null().sum())
        s = s.drop_nulls().cast(pl.Float64)
        if len(s) == 0:
            return {"mean": None, "sd": None, "n_missing": n_missing}
        return {
            "mean": round(float(s.mean()), 1),
            "sd": round(float(s.std()), 1),
            "n_missing": n_missing,
        }

    def _n_pct_raw(mask, n_total):
        n = int(mask.drop_nulls().sum())
        return {"n": n, "pct": round(n / n_total * 100, 1) if n_total else None, "N": n_total}

    def _cat_rows_raw(s, n_total):
        s = s.drop_nulls()
        if len(s) == 0:
            return {}
        vc = s.value_counts().sort("count", descending=True)
        col = [c for c in vc.columns if c != "count"][0]
        return {
            str(r[col]): {"n": int(r["count"]), "pct": round(r["count"] / n_total * 100, 1)}
            for r in vc.iter_rows(named=True)
        }

    def _table1_col(sub):
        n = sub.height
        sub_pt = sub.unique(subset=["hospitalization_id"])
        n_pt = sub_pt.height
        rows = [("N (vent-days)", f"{n:,}")]
        rows.append(("N (patients)", f"{n_pt:,}"))
        rows.append(("Age, years", _median_iqr(sub_pt["age_at_admission"])))
        rows.append((
            "Sex, female",
            _n_pct(sub_pt["sex_category"].str.to_lowercase() == "female", n_pt),
        ))
        rows.append(("Race", ""))
        rows.extend(_cat_rows(sub_pt["race_category"], n_pt))
        rows.append((
            "Ethnicity, Hispanic/Latino",
            _n_pct(
                sub_pt["ethnicity_category"].str.to_lowercase().str.starts_with("hispanic"),
                n_pt,
            ),
        ))
        _lang = sub_pt["language_category"].map_elements(
            lambda x: "English" if x and x.lower() == "english" else "Non-English",
            return_dtype=pl.Utf8,
        )
        rows.append(("Preferred language, English", _n_pct(_lang == "English", n_pt)))
        rows.append(("Weight, kg", _median_iqr(sub_pt["weight_kg"])))
        rows.append(("Weight missing, n", _n_pct(sub_pt["weight_kg"].is_null(), n_pt)))
        rows.append(("Height, cm", _median_iqr(sub_pt["height_cm"])))
        rows.append(("Height missing, n", _n_pct(sub_pt["height_cm"].is_null(), n_pt)))
        rows.append(("BMI, kg/m²", _median_iqr(sub_pt["bmi"])))
        rows.append(("BMI missing, n", _n_pct(sub_pt["bmi"].is_null(), n_pt)))
        rows.append(("CCI", _median_iqr(sub_pt["cci_score"])))
        rows.append(("SOFA prior day", _median_iqr(sub["sofa_prior"])))
        rows.append(("GCS prior day", _median_iqr(sub["gcs_prior"])))
        rows.append(("BZD exposure prior day", _n_pct(sub["bzd_prior"], n)))
        rows.append(("Propofol exposure prior day", _n_pct(sub["propofol_prior"], n)))
        rows.append(("Opioid infusion prior day", _n_pct(sub["opioid_prior"], n)))
        rows.append(("NMB prior day", _n_pct(sub["nmb_prior"], n)))
        rows.append(("RASS prior day", _median_iqr(sub["rass_prior"])))
        rows.append(("NEE prior day", _mean_sd(sub["nee_prior"])))
        rows.append(("FiO2 prior day", _mean_sd(sub["fio2_prior"])))
        rows.append(("PEEP prior day", _mean_sd(sub["peep_prior"])))
        rows.append(("Hospital LOS, days", _median_iqr(sub_pt["hospital_los_hours"] / 24)))
        rows.append(("ICU LOS, days", _median_iqr(sub_pt["first_icu_los_hours"] / 24)))
        rows.append(("IMV duration, hours", _median_iqr(sub_pt["imv_duration_hours"])))
        rows.append(("In-hospital mortality", _n_pct(
            sub_pt["discharge_category"].str.to_lowercase() == "expired", n_pt
        )))
        return rows

    def _table1_col_raw(sub):
        n = sub.height
        sub_pt = sub.unique(subset=["hospitalization_id"])
        n_pt = sub_pt.height
        d = {}
        d["N (vent-days)"] = {"n": n}
        d["N (patients)"] = {"n": n_pt}
        d["Age, years"] = _median_iqr_raw(sub_pt["age_at_admission"])
        d["Sex, female"] = _n_pct_raw(
            sub_pt["sex_category"].str.to_lowercase() == "female", n_pt
        )
        d["Race"] = _cat_rows_raw(sub_pt["race_category"], n_pt)
        d["Ethnicity, Hispanic/Latino"] = _n_pct_raw(
            sub_pt["ethnicity_category"].str.to_lowercase().str.starts_with("hispanic"), n_pt
        )
        _lang = sub_pt["language_category"].map_elements(
            lambda x: "English" if x and x.lower() == "english" else "Non-English",
            return_dtype=pl.Utf8,
        )
        d["Preferred language, English"] = _n_pct_raw(_lang == "English", n_pt)
        d["Weight, kg"] = _median_iqr_raw(sub_pt["weight_kg"])
        d["Weight missing, n"] = _n_pct_raw(sub_pt["weight_kg"].is_null(), n_pt)
        d["Height, cm"] = _median_iqr_raw(sub_pt["height_cm"])
        d["Height missing, n"] = _n_pct_raw(sub_pt["height_cm"].is_null(), n_pt)
        d["BMI, kg/m²"] = _median_iqr_raw(sub_pt["bmi"])
        d["BMI missing, n"] = _n_pct_raw(sub_pt["bmi"].is_null(), n_pt)
        d["CCI"] = _median_iqr_raw(sub_pt["cci_score"])
        d["SOFA prior day"] = _median_iqr_raw(sub["sofa_prior"])
        d["GCS prior day"] = _median_iqr_raw(sub["gcs_prior"])
        d["BZD exposure prior day"] = _n_pct_raw(sub["bzd_prior"], n)
        d["Propofol exposure prior day"] = _n_pct_raw(sub["propofol_prior"], n)
        d["Opioid infusion prior day"] = _n_pct_raw(sub["opioid_prior"], n)
        d["NMB prior day"] = _n_pct_raw(sub["nmb_prior"], n)
        d["RASS prior day"] = _median_iqr_raw(sub["rass_prior"])
        d["NEE prior day"] = _mean_sd_raw(sub["nee_prior"])
        d["FiO2 prior day"] = _mean_sd_raw(sub["fio2_prior"])
        d["PEEP prior day"] = _mean_sd_raw(sub["peep_prior"])
        d["Hospital LOS, days"] = _median_iqr_raw(sub_pt["hospital_los_hours"] / 24)
        d["ICU LOS, days"] = _median_iqr_raw(sub_pt["first_icu_los_hours"] / 24)
        d["IMV duration, hours"] = _median_iqr_raw(sub_pt["imv_duration_hours"])
        d["In-hospital mortality"] = _n_pct_raw(
            sub_pt["discharge_category"].str.to_lowercase() == "expired", n_pt
        )
        return d

    # ── Step 7: Compute for 3 groups and save ──
    _ov = _table1_col(_tbl)
    _el = dict(_table1_col(_tbl.filter(pl.col("eligible_event") == 1)))
    _inel = dict(_table1_col(_tbl.filter(
        pl.col("eligible_event").is_null() | (pl.col("eligible_event") != 1)
    )))

    _chars = [r[0] for r in _ov]
    table1 = pd.DataFrame({
        "Characteristic": _chars,
        "Overall": [r[1] for r in _ov],
        "SBT Eligible": [_el.get(c, "—") for c in _chars],
        "SBT Ineligible": [_inel.get(c, "—") for c in _chars],
    })

    table1.to_csv(OUTPUT_SHARE / "table1.csv", index=False)

    # ── Step 8: Export JSON for cross-site aggregation ──
    table1_dict = {
        "site": SITE_NAME,
        "groups": {
            "Overall": _table1_col_raw(_tbl),
            "SBT Eligible": _table1_col_raw(
                _tbl.filter(pl.col("eligible_event") == 1)
            ),
            "SBT Ineligible": _table1_col_raw(
                _tbl.filter(
                    pl.col("eligible_event").is_null()
                    | (pl.col("eligible_event") != 1)
                )
            ),
        },
    }
    with open(OUTPUT_SHARE / "table1.json", "w") as f:
        json.dump(table1_dict, f, indent=2)

    print(table1.to_string(index=False))
    return


@app.cell
def _(
    OUTPUT_SHARE,
    SITE_NAME,
    cohort,
    cohort_with_delivery,
    df_grouped,
    n_vent_days,
    pl,
    sbt_flag_cols,
):
    _all_flags = sbt_flag_cols + ["sbt_ground_truth"]

    if "hospital_id" in cohort.columns:
        # Multi-hospital: join hospital_id onto df_grouped, then summarize per hospital
        _hosp_ids = cohort.select("hospitalization_id", "hospital_id").unique()
        _hid = cohort_with_delivery.select("hospitalization_id", "hosp_id_icu_day").unique()
        _final = df_grouped.join(_hid, on="hosp_id_icu_day", how="left")
        _final = _final.join(_hosp_ids, on="hospitalization_id", how="left")

        _agg_exprs = [
            (pl.col("eligible_event") == 1).sum().cast(pl.Int64).alias("eligible_event_count"),
        ]
        for _flag in _all_flags:
            _agg_exprs.append(
                (pl.col(_flag) == 1).sum().cast(pl.Int64).alias(f"{_flag}_count")
            )

        hospital_summary = (
            _final.filter(pl.col("hospital_id").is_not_null())
            .group_by("hospital_id")
            .agg(_agg_exprs)
        )
        hospital_summary = hospital_summary.with_columns(
            pl.lit(n_vent_days).alias("total_vent_days"),
        )
        hospital_summary = hospital_summary.with_columns(
            pl.concat_str([pl.lit(SITE_NAME + "_"), pl.col("hospital_id").cast(pl.Utf8)]).alias("Site_Hospital"),
        )
        for _flag in _all_flags:
            hospital_summary = hospital_summary.with_columns(
                pl.when(pl.col("eligible_event_count") > 0)
                .then((pl.col(f"{_flag}_count") / pl.col("eligible_event_count") * 100).round(2))
                .otherwise(0.0)
                .alias(f"pct_{_flag}")
            )
    else:
        # Single-site: one summary row for the whole site
        _eligible_n = df_grouped.filter(pl.col("eligible_event") == 1).height
        _row = {"Site_Hospital": SITE_NAME, "total_vent_days": n_vent_days, "eligible_event_count": _eligible_n}
        for _flag in _all_flags:
            _n = df_grouped.filter(pl.col(_flag) == 1).height
            _row[f"{_flag}_count"] = _n
            _row[f"pct_{_flag}"] = round(_n / _eligible_n * 100, 2) if _eligible_n > 0 else 0.0
        hospital_summary = pl.DataFrame([_row])

    hospital_summary.write_csv(OUTPUT_SHARE / f"sbt_stats_{SITE_NAME}.csv")
    print(hospital_summary.to_pandas().T)
    return


@app.cell
def _(OUTPUT_SBT, df_grouped):
    _path = OUTPUT_SBT / "sbt_day_level.parquet"
    df_grouped.write_parquet(str(_path))
    print(f"Saved: {_path} ({df_grouped.height:,} rows x {df_grouped.width} cols)")
    return


@app.cell
def _(df_grouped):
    df_grouped['eligibility_failure_reason'].value_counts()
    return


@app.cell
def _(df_grouped):
    df_grouped['hosp_id_icu_day'].n_unique()
    return


@app.cell
def _(df):
    df
    return


if __name__ == "__main__":
    app.run()
