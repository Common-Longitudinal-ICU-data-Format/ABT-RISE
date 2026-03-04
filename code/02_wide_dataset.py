import marimo

__generated_with = "0.20.2"
app = marimo.App(width="full")


@app.cell
def _():
    import polars as pl
    import json
    from pathlib import Path

    _config_path = Path(__file__).parent.parent / "clif_config.json"
    with open(_config_path) as _f:
        _cfg = json.load(_f)

    OUTPUT_PHI = Path(__file__).parent.parent / "output_phi"

    ALL_MEDS = [
        "norepinephrine", "epinephrine", "phenylephrine", "angiotensin",
        "vasopressin", "dopamine", "dobutamine", "milrinone", "isoproterenol",
        "cisatracurium", "vecuronium", "rocuronium",
        "fentanyl", "propofol", "lorazepam", "midazolam",
        "hydromorphone", "morphine",
    ]

    ALL_ASSESSMENTS = [
        "sat_screen_pass_fail", "sat_screen_performed",
        "sat_delivery_pass_fail", "sat_delivery_performed",
        "sbt_screen_pass_fail", "sbt_screen_performed",
        "sbt_delivery_pass_fail", "sbt_delivery_performed",
        "rass", "gcs_total",
    ]

    print(f"Output directory: {OUTPUT_PHI}")
    print(f"Target meds: {len(ALL_MEDS)}, Target assessments: {len(ALL_ASSESSMENTS)}")
    return ALL_ASSESSMENTS, ALL_MEDS, OUTPUT_PHI, pl


@app.cell
def _(OUTPUT_PHI, pl):
    cohort = pl.read_parquet(
        OUTPUT_PHI / "cohort.parquet",
        columns=["hospitalization_id", "first_icu_start", "first_icu_end"],
    )
    cohort_ids = cohort.select("hospitalization_id")
    print(f"Cohort: {cohort.height:,} hospitalizations")
    print(f"ICU start range: {cohort['first_icu_start'].min()} → {cohort['first_icu_start'].max()}")
    return cohort, cohort_ids


@app.cell
def _(OUTPUT_PHI, cohort, cohort_ids, pl):
    _resp_cols = [
        "hospitalization_id", "recorded_dttm", "is_scaffold",
        "device_category", "mode_category", "device_name", "mode_name",
        "fio2_set", "peep_set", "resp_rate_set", "pressure_support_set",
    ]
    _resp_raw = pl.read_parquet(OUTPUT_PHI / "resp_waterfall_cohort.parquet", columns=_resp_cols)

    # Filter to cohort IDs
    _resp = _resp_raw.join(cohort_ids, on="hospitalization_id", how="semi")
    print(f"Resp records after cohort filter: {_resp.height:,} (from {_resp_raw.height:,})")

    # Strip timezone (wall-clock time matches naive cohort datetimes)
    _resp = _resp.with_columns(pl.col("recorded_dttm").dt.replace_time_zone(None))

    # Keep all rows (scaffold + real), drop the is_scaffold column
    _resp = _resp.drop("is_scaffold")

    # Filter to ICU window (exact datetimes, no truncation)
    _resp = (
        _resp
        .join(cohort.select("hospitalization_id", "first_icu_start", "first_icu_end"),
              on="hospitalization_id", how="inner")
        .filter(
            (pl.col("recorded_dttm") >= pl.col("first_icu_start"))
            & (pl.col("recorded_dttm") <= pl.col("first_icu_end"))
        )
        .drop(["first_icu_start", "first_icu_end"])
    )

    # Deduplicate per (hosp, recorded_dttm): last() for categoricals, mean() for numerics
    _cat_cols = ["device_category", "mode_category", "device_name", "mode_name"]
    _num_cols = ["fio2_set", "peep_set", "resp_rate_set", "pressure_support_set"]

    resp_df = (
        _resp
        .sort(["hospitalization_id", "recorded_dttm"])
        .group_by(["hospitalization_id", "recorded_dttm"])
        .agg(
            [pl.col(c).last() for c in _cat_cols]
            + [pl.col(c).mean() for c in _num_cols]
        )
    )
    print(f"Resp (original timestamps): {resp_df.height:,} rows")
    return (resp_df,)


@app.cell
def _(OUTPUT_PHI, cohort, cohort_ids, pl):
    _spo2_raw = pl.read_parquet(
        OUTPUT_PHI / "vitals_spo2_cohort.parquet",
        columns=["hospitalization_id", "recorded_dttm", "vital_value"],
    )

    # Filter to cohort IDs
    _spo2 = _spo2_raw.join(cohort_ids, on="hospitalization_id", how="semi")
    print(f"SpO2 records after cohort filter: {_spo2.height:,} (from {_spo2_raw.height:,})")

    # Filter to ICU window (exact datetimes, no truncation)
    _spo2 = (
        _spo2
        .join(cohort.select("hospitalization_id", "first_icu_start", "first_icu_end"),
              on="hospitalization_id", how="inner")
        .filter(
            (pl.col("recorded_dttm") >= pl.col("first_icu_start"))
            & (pl.col("recorded_dttm") <= pl.col("first_icu_end"))
        )
        .drop(["first_icu_start", "first_icu_end"])
    )

    # Deduplicate per (hosp, recorded_dttm) via mean — no truncation
    spo2_df = (
        _spo2
        .group_by(["hospitalization_id", "recorded_dttm"])
        .agg(pl.col("vital_value").mean().alias("spo2"))
    )
    print(f"SpO2 (original timestamps): {spo2_df.height:,} rows")
    return (spo2_df,)


@app.cell
def _(ALL_MEDS, OUTPUT_PHI, cohort, pl):
    _meds_raw = pl.read_parquet(
        OUTPUT_PHI / "meds_cohort.parquet",
        columns=["hospitalization_id", "admin_dttm", "med_category", "med_dose_converted"],
    )
    print(f"Meds records: {_meds_raw.height:,}")

    # Filter to ICU window (exact datetimes, no truncation)
    _meds = (
        _meds_raw
        .join(cohort.select("hospitalization_id", "first_icu_start", "first_icu_end"),
              on="hospitalization_id", how="inner")
        .filter(
            (pl.col("admin_dttm") >= pl.col("first_icu_start"))
            & (pl.col("admin_dttm") <= pl.col("first_icu_end"))
        )
        .drop(["first_icu_start", "first_icu_end"])
    )

    # Rename admin_dttm → recorded_dttm (no truncation)
    _meds = _meds.rename({"admin_dttm": "recorded_dttm"})

    # Aggregate: mean dose per (hosp, recorded_dttm, med_category)
    _meds_agg = (
        _meds
        .group_by(["hospitalization_id", "recorded_dttm", "med_category"])
        .agg(pl.col("med_dose_converted").mean())
    )

    # Pivot to wide
    meds_wide = _meds_agg.pivot(
        on="med_category",
        index=["hospitalization_id", "recorded_dttm"],
        values="med_dose_converted",
    )

    # Add null columns for any missing med categories
    _existing = set(meds_wide.columns) - {"hospitalization_id", "recorded_dttm"}
    for med in ALL_MEDS:
        if med not in _existing:
            meds_wide = meds_wide.with_columns(pl.lit(None).cast(pl.Float64).alias(med))

    # Reorder to canonical order
    meds_wide = meds_wide.select(["hospitalization_id", "recorded_dttm"] + ALL_MEDS)

    # -------------------------------------------------------------------------
    # Compute Norepinephrine Equivalent (NEE) from vasopressors
    # -------------------------------------------------------------------------
    _nee_weights = {
        "norepinephrine": 1.0,
        "epinephrine": 1.0,
        "phenylephrine": 1 / 10.0,
        "dopamine": 1 / 100.0,
        "vasopressin": 2.5,
        "angiotensin": 10.0,
    }
    _vaso_cols = [c for c in _nee_weights if c in meds_wide.columns]
    _nee_expr = pl.sum_horizontal(
        [pl.col(c).fill_null(0) * w for c, w in _nee_weights.items() if c in meds_wide.columns]
    )
    # NEE is null when ALL vasopressor columns are null (no vasopressors recorded)
    _all_vaso_null = pl.all_horizontal([pl.col(c).is_null() for c in _vaso_cols])
    meds_wide = meds_wide.with_columns(
        pl.when(_all_vaso_null).then(None).otherwise(_nee_expr).cast(pl.Float64).alias("nee")
    )

    print(f"Meds wide (original timestamps): {meds_wide.height:,} rows x {meds_wide.width} cols")
    print(f"  Missing categories (added as null): {sorted(_existing.symmetric_difference(set(ALL_MEDS)) - _existing)}")
    print(f"  NEE computed: {meds_wide['nee'].is_not_null().sum():,} non-null values")

    # -------------------------------------------------------------------------
    # Consolidated sedation & paralytic dose columns (SAT notebook approach)
    # -------------------------------------------------------------------------
    _sedatives = ["propofol", "lorazepam", "midazolam"]
    _opioids = ["fentanyl", "morphine", "hydromorphone"]
    _all_sedation = _sedatives + _opioids
    _paralytics = ["cisatracurium", "vecuronium", "rocuronium"]

    _sed_cols = [c for c in _all_sedation if c in meds_wide.columns]
    _sed_only_cols = [c for c in _sedatives if c in meds_wide.columns]
    _par_cols = [c for c in _paralytics if c in meds_wide.columns]

    meds_wide = meds_wide.with_columns(
        # Min dose across all 6 sedation meds (nulls skipped)
        pl.min_horizontal([pl.col(c) for c in _sed_cols]).alias("min_sedation_dose"),
        # Min dose across all 6 where dose > 0 (zeros treated as null)
        pl.min_horizontal(
            [pl.when(pl.col(c) > 0).then(pl.col(c)) for c in _sed_cols]
        ).alias("min_sedation_dose_2"),
        # Min dose across sedatives only (no opioids), null → 0
        pl.min_horizontal([pl.col(c) for c in _sed_only_cols])
        .fill_null(0)
        .alias("min_sedation_dose_non_ops"),
        # Max dose across paralytics, null → 0
        pl.max_horizontal([pl.col(c) for c in _par_cols])
        .fill_null(0)
        .alias("max_paralytics"),
    )
    print(
        f"  Sedation consolidated: "
        f"min_sedation_dose={meds_wide['min_sedation_dose'].is_not_null().sum():,}, "
        f"max_paralytics={( meds_wide['max_paralytics'] > 0 ).sum():,} non-zero"
    )
    return (meds_wide,)


@app.cell
def _(ALL_ASSESSMENTS, OUTPUT_PHI, cohort, pl):
    _assess_raw = pl.read_parquet(
        OUTPUT_PHI / "assessments_cohort.parquet",
        columns=["hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value"],
    )
    print(f"Assessment records: {_assess_raw.height:,}")

    # Filter to ICU window (exact datetimes, no truncation)
    _assess = (
        _assess_raw
        .join(cohort.select("hospitalization_id", "first_icu_start", "first_icu_end"),
              on="hospitalization_id", how="inner")
        .filter(
            (pl.col("recorded_dttm") >= pl.col("first_icu_start"))
            & (pl.col("recorded_dttm") <= pl.col("first_icu_end"))
        )
        .drop(["first_icu_start", "first_icu_end"])
    )

    # Lowercase assessment_category
    _assess = _assess.with_columns(pl.col("assessment_category").str.to_lowercase())

    # Aggregate: mean per (hosp, recorded_dttm, assessment_category) — no truncation
    _assess_agg = (
        _assess
        .group_by(["hospitalization_id", "recorded_dttm", "assessment_category"])
        .agg(pl.col("numerical_value").mean())
    )

    # Pivot to wide
    assess_wide = _assess_agg.pivot(
        on="assessment_category",
        index=["hospitalization_id", "recorded_dttm"],
        values="numerical_value",
    )

    # Add null columns for any missing assessment categories
    _existing = set(assess_wide.columns) - {"hospitalization_id", "recorded_dttm"}
    for cat in ALL_ASSESSMENTS:
        if cat not in _existing:
            assess_wide = assess_wide.with_columns(pl.lit(None).cast(pl.Float64).alias(cat))

    # Reorder to canonical order
    assess_wide = assess_wide.select(["hospitalization_id", "recorded_dttm"] + ALL_ASSESSMENTS)
    print(f"Assessments wide (original timestamps): {assess_wide.height:,} rows x {assess_wide.width} cols")
    print(f"  Missing categories (added as null): {sorted(set(ALL_ASSESSMENTS) - _existing)}")
    return (assess_wide,)


@app.cell
def _(assess_wide, cohort, meds_wide, pl, resp_df, spo2_df):
    # Build spine: union of all unique (hospitalization_id, recorded_dttm) from 4 sources
    _key_cols = ["hospitalization_id", "recorded_dttm"]
    spine = (
        pl.concat([
            resp_df.select(_key_cols),
            spo2_df.select(_key_cols),
            meds_wide.select(_key_cols),
            assess_wide.select(_key_cols),
        ])
        .unique()
    )
    print(f"Spine (union of all timestamps): {spine.height:,} rows")

    # Join cohort metadata onto spine
    spine = spine.join(
        cohort.select("hospitalization_id", "first_icu_start", "first_icu_end"),
        on="hospitalization_id",
        how="left",
    )

    # Left join each source onto the spine
    wide = (
        spine
        .join(resp_df, on=_key_cols, how="left")
        .join(spo2_df, on=_key_cols, how="left")
        .join(meds_wide, on=_key_cols, how="left")
        .join(assess_wide, on=_key_cols, how="left")
        .sort(_key_cols)
    )

    # Calendar-based ICU day (day 1 = date of first_icu_start)
    wide = wide.with_columns(
        ((pl.col("recorded_dttm").dt.date() - pl.col("first_icu_start").dt.date()).dt.total_days() + 1)
        .cast(pl.Int32)
        .alias("icu_day")
    )

    # Vent-day columns (shared by SAT & SBT notebooks)
    VENT_DAY_ANCHOR_HOUR = 0  # calendar day (midnight-to-midnight)

    wide = wide.with_columns(
        (pl.col("recorded_dttm") - pl.duration(hours=VENT_DAY_ANCHOR_HOUR))
        .dt.truncate("1d")
        .alias("vent_day_date")
    )

    wide = wide.with_columns(
        pl.col("vent_day_date")
        .rank("dense")
        .over("hospitalization_id")
        .cast(pl.Int64)
        .alias("vent_day_num")
    )

    wide = wide.with_columns(
        pl.concat_str([
            pl.col("hospitalization_id").cast(pl.Utf8),
            pl.lit("_day_"),
            pl.col("vent_day_num").cast(pl.Utf8),
        ]).alias("hosp_id_day_key")
    )

    print(f"Wide dataset: {wide.height:,} rows x {wide.width} cols")
    print(f"\nNull counts (before forward-fill):")
    _nulls = wide.null_count()
    for col in wide.columns:
        _n = _nulls[col][0]
        _pct = _n / wide.height * 100
        print(f"  {col:30s} {_n:>10,} ({_pct:5.1f}%)")

    # Forward-fill within each patient (all columns except device_name, mode_name)
    _no_fill = {
        "hospitalization_id", "recorded_dttm", "device_name", "mode_name", "icu_day",
        "vent_day_date", "vent_day_num", "hosp_id_day_key",
        # SAT/SBT flowsheet flags are point-in-time; forward-fill bleeds across days
        "sat_screen_pass_fail", "sat_screen_performed",
        "sat_delivery_pass_fail", "sat_delivery_performed",
        "sbt_screen_pass_fail", "sbt_screen_performed",
        "sbt_delivery_pass_fail", "sbt_delivery_performed",
    }
    _fill_cols = [c for c in wide.columns if c not in _no_fill]

    wide = wide.with_columns(
        pl.col(_fill_cols).forward_fill().over("hospitalization_id")
    )

    print(f"\nNull counts (after forward-fill):")
    _nulls = wide.null_count()
    for col in wide.columns:
        _n = _nulls[col][0]
        _pct = _n / wide.height * 100
        print(f"  {col:30s} {_n:>10,} ({_pct:5.1f}%)")
    return (wide,)


@app.cell
def _(OUTPUT_PHI, wide):
    _path = OUTPUT_PHI / "wide_dataset.parquet"
    wide.write_parquet(str(_path))
    _size_mb = _path.stat().st_size / 1024 / 1024
    print(f"Saved: {_path}")
    print(f"Shape: {wide.height:,} rows x {wide.width} cols")
    print(f"File size: {_size_mb:.1f} MB")
    return


@app.cell
def _(wide):
    wide
    return


if __name__ == "__main__":
    app.run()
