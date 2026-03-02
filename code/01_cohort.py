import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import polars as pl
    import pandas as pd
    import json
    from datetime import timedelta
    from pathlib import Path

    def strip_tz(df: pd.DataFrame) -> pd.DataFrame:
        """Remove timezone from all tz-aware datetime columns, preserving wall-clock times."""
        for col in df.select_dtypes(include=["datetimetz"]).columns:
            df[col] = df[col].dt.tz_localize(None)
        return df

    _config_path = Path(__file__).parent.parent / "clif_config.json"
    with open(_config_path) as _f:
        _cfg = json.load(_f)

    SITE_NAME = _cfg["site_name"]
    DATA_DIR = _cfg["data_directory"]
    FILETYPE = _cfg["filetype"]
    TIMEZONE = _cfg["timezone"]
    OUTPUT_DIR = Path(__file__).parent.parent / "output_to_share"
    OUTPUT_PHI = Path(__file__).parent.parent / "output_phi"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_PHI.mkdir(parents=True, exist_ok=True)
    return (
        DATA_DIR,
        FILETYPE,
        OUTPUT_DIR,
        OUTPUT_PHI,
        SITE_NAME,
        TIMEZONE,
        json,
        pd,
        pl,
        strip_tz,
        timedelta,
    )


@app.cell
def _(DATA_DIR, FILETYPE, TIMEZONE, pl, strip_tz):
    from clifpy import Hospitalization

    _hosp = Hospitalization.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE, verbose=False,
    )
    hosp_df = pl.from_pandas(strip_tz(_hosp.df))
    n_total = len(hosp_df)
    print(f"Total hospitalizations: {n_total:,}")
    return hosp_df, n_total


@app.cell
def _(DATA_DIR, FILETYPE, TIMEZONE, pl, strip_tz):
    from clifpy import Patient

    _patient = Patient.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE, verbose=False,
    )
    patient_df = pl.from_pandas(strip_tz(_patient.df))
    print(f"Total patients: {len(patient_df):,}")
    return (patient_df,)


@app.cell
def _(hosp_df, pl):
    hosp_adults = hosp_df.filter(pl.col("age_at_admission") >= 18)
    n_excluded_age = len(hosp_df) - len(hosp_adults)
    print(f"Excluded (age < 18): {n_excluded_age:,}")
    print(f"Remaining: {len(hosp_adults):,}")
    return hosp_adults, n_excluded_age


@app.cell
def _(SITE_NAME, hosp_adults, pl):
    if SITE_NAME.upper() == "MIMIC":
        hosp_dated = hosp_adults
        n_excluded_dates = 0
        print("MIMIC site detected — skipping 2018-2024 date filter")
    else:
        hosp_dated = hosp_adults.filter(
            (pl.col("admission_dttm").dt.year() >= 2018)
            & (pl.col("admission_dttm").dt.year() <= 2024)
            & (pl.col("discharge_dttm").dt.year() <= 2024)
        )
        n_excluded_dates = len(hosp_adults) - len(hosp_dated)
    print(f"Excluded (outside 2018-2024): {n_excluded_dates:,}")
    print(f"Remaining: {len(hosp_dated):,}")
    return hosp_dated, n_excluded_dates


@app.cell
def _(DATA_DIR, FILETYPE, TIMEZONE, hosp_dated, pl, strip_tz):
    from clifpy import Adt

    _adt = Adt.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE, verbose=False,
    )
    adt_df = pl.from_pandas(strip_tz(_adt.df))

    _icu_hosp_ids = (
        adt_df
        .filter(pl.col("location_category").str.to_lowercase() == "icu")
        .select("hospitalization_id")
        .unique()
    )
    hosp_icu = hosp_dated.join(_icu_hosp_ids, on="hospitalization_id", how="inner")
    n_excluded_no_icu = len(hosp_dated) - len(hosp_icu)
    print(f"Excluded (no ICU stay): {n_excluded_no_icu:,}")
    print(f"Remaining: {len(hosp_icu):,}")
    return adt_df, hosp_icu, n_excluded_no_icu


@app.cell
def _(adt_df, hosp_icu, pl):
    _cohort_ids = hosp_icu.select("hospitalization_id")
    _adt_cohort = adt_df.join(_cohort_ids, on="hospitalization_id", how="inner")
    _adt_sorted = _adt_cohort.sort(["hospitalization_id", "in_dttm"])

    # Flag ICU/procedural rows, detect transitions, assign block IDs
    _adt_flagged = (
        _adt_sorted
        .with_columns(
            pl.col("location_category")
            .str.to_lowercase().is_in(["icu", "procedural"])
            .alias("is_icu_proc")
        )
        .with_columns(
            (
                pl.col("is_icu_proc")
                != pl.col("is_icu_proc").shift(1).over("hospitalization_id")
            )
            .fill_null(True)
            .alias("transition")
        )
        .with_columns(
            pl.col("transition")
            .cum_sum()
            .over("hospitalization_id")
            .alias("block_id")
        )
    )

    # Merge each contiguous ICU/procedural block into one stay
    _icu_blocks = (
        _adt_flagged
        .filter(pl.col("is_icu_proc"))
        .group_by(["hospitalization_id", "block_id"])
        .agg(
            pl.col("in_dttm").min().alias("icu_in_dttm"),
            pl.col("out_dttm").max().alias("icu_out_dttm"),
            (pl.col("location_category") == "icu").any().alias("has_icu"),
        )
        .filter(pl.col("has_icu"))
        .drop(["block_id", "has_icu"])
    )

    # Drop stays with null icu_out_dttm
    icu_stays = (
        _icu_blocks
        .filter(pl.col("icu_out_dttm").is_not_null())
        .sort(["hospitalization_id", "icu_in_dttm"])
        .with_columns(
            pl.col("icu_in_dttm")
            .rank("ordinal")
            .over("hospitalization_id")
            .cast(pl.UInt32)
            .alias("icu_rank")
        )
    )
    _valid_hosp_ids = icu_stays.select("hospitalization_id").unique()
    hosp_after_merge = hosp_icu.join(
        _valid_hosp_ids, on="hospitalization_id", how="inner"
    )
    n_excluded_null_out = len(hosp_icu) - len(hosp_after_merge)
    print(f"Excluded (all ICU stays had null out_dttm): {n_excluded_null_out:,}")
    print(f"Remaining: {len(hosp_after_merge):,}")
    print(f"Total merged ICU stays: {len(icu_stays):,}")
    return hosp_after_merge, icu_stays, n_excluded_null_out


@app.cell
def _(hosp_after_merge, icu_stays, pl, timedelta):
    # Keep ALL ICU stays >= 24h (not just the first)
    icu_24h = (
        icu_stays
        .with_columns(
            (pl.col("icu_out_dttm") - pl.col("icu_in_dttm")).alias("icu_los")
        )
        .filter(pl.col("icu_los") >= timedelta(hours=24))
    )
    _valid = icu_24h.select("hospitalization_id").unique()
    hosp_24h = hosp_after_merge.join(_valid, on="hospitalization_id", how="inner")
    n_excluded_short_icu = len(hosp_after_merge) - len(hosp_24h)
    print(f"Excluded (no ICU stay >= 24h): {n_excluded_short_icu:,}")
    print(f"Remaining: {len(hosp_24h):,}")
    return hosp_24h, icu_24h, n_excluded_short_icu


@app.cell
def _(DATA_DIR, FILETYPE, TIMEZONE, icu_24h, pl, strip_tz):
    from clifpy import RespiratorySupport

    _rs = RespiratorySupport.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE, verbose=False,
    )
    _rs_all = pl.from_pandas(strip_tz(_rs.df))

    # Join with ALL qualifying ICU stays (>= 24h), not just the first
    rs_df = (
        _rs_all
        .join(
            icu_24h.select(["hospitalization_id", "icu_in_dttm", "icu_out_dttm"]),
            on="hospitalization_id",
            how="inner",
        )
        .filter(
            (pl.col("recorded_dttm") >= pl.col("icu_in_dttm"))
            & (pl.col("recorded_dttm") <= pl.col("icu_out_dttm"))
        )
    )
    print(f"Respiratory support records in qualifying ICU windows: {len(rs_df):,}")
    return (rs_df,)


@app.cell
def _(hosp_24h, icu_24h, pl, rs_df):
    # --- Exclusion: tracheostomy == 1 at any point in ICU windows ---
    _trach_ids = (
        rs_df.filter(pl.col("tracheostomy") == 1)
        .select("hospitalization_id")
        .unique()
    )
    _hosp_no_trach = hosp_24h.join(_trach_ids, on="hospitalization_id", how="anti")
    _icu_no_trach = icu_24h.join(_trach_ids, on="hospitalization_id", how="anti")
    n_excluded_trach = len(hosp_24h) - len(_hosp_no_trach)
    print(f"Excluded (tracheostomy): {n_excluded_trach:,}")
    print(f"Remaining: {len(_hosp_no_trach):,}")

    # --- Exclusion: Trach Collar device in ICU windows ---
    _collar_ids = (
        rs_df.filter(pl.col("device_category").str.to_lowercase() == "trach collar")
        .select("hospitalization_id")
        .unique()
    )
    _hosp_no_collar = _hosp_no_trach.join(_collar_ids, on="hospitalization_id", how="anti")
    _icu_no_collar = _icu_no_trach.join(_collar_ids, on="hospitalization_id", how="anti")
    n_excluded_collar = len(_hosp_no_trach) - len(_hosp_no_collar)
    print(f"Excluded (trach collar): {n_excluded_collar:,}")
    print(f"Remaining: {len(_hosp_no_collar):,}")

    # --- Find first ICU stay with IMV (from remaining) ---
    _imv_stays = (
        rs_df.filter(pl.col("device_category").str.to_lowercase() == "imv")
        .select(["hospitalization_id", "icu_in_dttm"])
        .unique()
    )
    first_icu = (
        _icu_no_collar
        .join(_imv_stays, on=["hospitalization_id", "icu_in_dttm"], how="semi")
        .sort(["hospitalization_id", "icu_in_dttm"])
        .unique(subset=["hospitalization_id"], keep="first")
    )
    hosp_final = _hosp_no_collar.join(
        first_icu.select("hospitalization_id"), on="hospitalization_id", how="inner"
    )
    n_excluded_no_imv = len(_hosp_no_collar) - len(hosp_final)
    print(f"Excluded (no IMV in any ICU stay >= 24h): {n_excluded_no_imv:,}")
    print(f"Final cohort: {len(hosp_final):,}")
    return (
        first_icu,
        hosp_final,
        n_excluded_collar,
        n_excluded_no_imv,
        n_excluded_trach,
    )


@app.cell
def _(DATA_DIR, FILETYPE, TIMEZONE, first_icu, hosp_final, pd, pl, strip_tz):
    import clifpy as _clifpy
    import pathlib as _pathlib

    _cohort_ids = hosp_final["hospitalization_id"].to_list()
    _cache_path = _pathlib.Path(__file__).parent.parent / "output_phi" / "resp_waterfall_cohort.parquet"

    if _cache_path.exists():
        print(f"Loading cached waterfall from: {_cache_path}")
        _resp_pd = pd.read_parquet(_cache_path)
        _resp_pd = _resp_pd[_resp_pd["hospitalization_id"].isin(_cohort_ids)].copy()
    else:
        _resp = _clifpy.RespiratorySupport.from_file(
            data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE,
            filters={"hospitalization_id": _cohort_ids}, verbose=False,
        )
        _resp.df["device_category"] = _resp.df["device_category"].str.lower()
        _resp_filled = _resp.waterfall(bfill=False, verbose=True)
        _resp_pd = _resp_filled.df.copy()
        _cache_path.parent.mkdir(parents=True, exist_ok=True)
        _resp_pd.to_parquet(_cache_path, index=False)
        print(f"Saved waterfall cache to: {_cache_path}")

    _resp_pd = strip_tz(_resp_pd)
    _resp_pl = pl.from_pandas(_resp_pd)

    # Window: admission_dttm → icu_out_dttm (first ICU with IMV)
    resp_wf = (
        _resp_pl
        .join(
            hosp_final.select(["hospitalization_id", "admission_dttm"]),
            on="hospitalization_id", how="inner",
        )
        .join(
            first_icu.select(["hospitalization_id", "icu_out_dttm"]),
            on="hospitalization_id", how="inner",
        )
        .filter(
            (pl.col("recorded_dttm") >= pl.col("admission_dttm"))
            & (pl.col("recorded_dttm") <= pl.col("icu_out_dttm"))
        )
        .sort(["hospitalization_id", "recorded_dttm"])
    )
    print(f"Waterfall resp records (admission→ICU end): {len(resp_wf):,}")
    return (resp_wf,)


@app.cell
def _(hosp_final, pl, resp_wf):
    # Forward-fill device, compute is_imv + lag/lead
    _rs_filled = (
        resp_wf
        .with_columns(
            pl.col("device_category")
            # .forward_fill()
            # .over("hospitalization_id")
            .alias("device_filled")
        )
        .with_columns(
            (pl.col("device_filled") == "imv").cast(pl.Int8).alias("is_imv")
        )
        .with_columns([
            pl.col("is_imv").shift(1).over("hospitalization_id").alias("lag1"),
            pl.col("is_imv").shift(2).over("hospitalization_id").alias("lag2"),
            pl.col("is_imv").shift(-1).over("hospitalization_id").alias("lead1"),
            pl.col("is_imv").shift(-2).over("hospitalization_id").alias("lead2"),
            pl.col("device_filled").shift(1).over("hospitalization_id").alias("lag1_device"),
            pl.col("device_filled").shift(-1).over("hospitalization_id").alias("lead1_device"),
        ])
        .with_columns(
            ((pl.col("is_imv") == 1) & pl.col("lpm_set").is_not_null() & (pl.col("lpm_set") > 0))
            .alias("imv_lpm_flag")
        )
        .with_columns(
            pl.col("imv_lpm_flag").shift(1).over("hospitalization_id").alias("lag1_imv_lpm")
        )
    )

    # --- Intubation: is_imv==1, (lag1==0|null), (lag2==0|null) ---
    _intubations = _rs_filled.filter(
        (pl.col("is_imv") == 1)
        & ((pl.col("lag1") == 0) | pl.col("lag1").is_null())
        & ((pl.col("lag2") == 0) | pl.col("lag2").is_null())
    ).sort(["hospitalization_id", "recorded_dttm"])

    _first_intubation = (
        _intubations
        .unique(subset=["hospitalization_id"], keep="first")
        .select([
            "hospitalization_id",
            pl.col("recorded_dttm").alias("intubation_time"),
            pl.col("lag1_device").alias("device_before_intubation"),
        ])
    )

    # --- Extubation: is_imv==1, lead1==0, lead2==0 (strict 2-row) ---
    _extubations = _rs_filled.filter(
        (pl.col("is_imv") == 1)
        & (pl.col("lead1") == 0)
        & (pl.col("lead2") == 0)
    ).sort(["hospitalization_id", "recorded_dttm"])

    _first_extubation = (
        _extubations
        .join(_first_intubation.select(["hospitalization_id", "intubation_time"]),
              on="hospitalization_id", how="inner")
        .filter(pl.col("recorded_dttm") >= pl.col("intubation_time"))
        .unique(subset=["hospitalization_id"], keep="first")
        .select([
            "hospitalization_id",
            pl.col("recorded_dttm").alias("extubation_time"),
            pl.col("lead1_device").alias("device_after_extubation"),
        ])
    )

    # --- Extubation via IMV+LPM co-charting (definitive device switch) ---
    _extub_lpm = (
        _rs_filled.filter(
            pl.col("imv_lpm_flag")
            & ((pl.col("lag1_imv_lpm") == False) | pl.col("lag1_imv_lpm").is_null())
        )
        .join(_first_intubation.select(["hospitalization_id", "intubation_time"]),
              on="hospitalization_id", how="inner")
        .filter(pl.col("recorded_dttm") >= pl.col("intubation_time"))
        .sort(["hospitalization_id", "recorded_dttm"])
        .unique(subset=["hospitalization_id"], keep="first")
        .select([
            "hospitalization_id",
            pl.col("recorded_dttm").alias("extubation_time"),
            pl.lit("Cannula or Facemask").alias("device_after_extubation"),
            pl.lit("2-row IMV+LPM").alias("extubation_method"),
        ])
    )

    # --- Combine both extubation methods (earliest per hospitalization) ---
    _first_extubation_combined = (
        pl.concat([
            _first_extubation.with_columns(pl.lit("2-row IMV").alias("extubation_method")),
            _extub_lpm,
        ])
        .sort(["hospitalization_id", "extubation_time"])
        .unique(subset=["hospitalization_id"], keep="first")
    )

    # --- No-extubation reason ---
    _hosp_no_extub = _first_intubation.join(
        _first_extubation_combined.select("hospitalization_id"),
        on="hospitalization_id", how="anti",
    ).select("hospitalization_id")

    _last_rs = (
        _rs_filled
        .join(_hosp_no_extub, on="hospitalization_id", how="semi")
        .unique(subset=["hospitalization_id"], keep="last")
        .select(["hospitalization_id", "is_imv"])
    )

    _unconfirmed = (
        _rs_filled
        .join(_hosp_no_extub, on="hospitalization_id", how="semi")
        .filter(
            (pl.col("is_imv") == 1) & (pl.col("lead1") == 0) & pl.col("lead2").is_null()
        )
        .select("hospitalization_id").unique()
    )

    _no_extub_classified = (
        _last_rs
        .join(hosp_final.select(["hospitalization_id", "discharge_category"]),
              on="hospitalization_id", how="left")
        .join(_unconfirmed.with_columns(pl.lit(True).alias("has_unconfirmed")),
              on="hospitalization_id", how="left")
        .with_columns(
            pl.when((pl.col("discharge_category").str.to_lowercase() == "expired") & (pl.col("is_imv") == 1))
            .then(pl.lit("Died on IMV"))
            .when(pl.col("has_unconfirmed") == True)
            .then(pl.lit("Unconfirmed extubation (single non-IMV reading)"))
            .when(pl.col("is_imv") == 1)
            .then(pl.lit("Still on IMV at ICU end"))
            .otherwise(pl.lit("Unknown"))
            .alias("no_extubation_reason")
        )
        .select(["hospitalization_id", "no_extubation_reason"])
    )

    # --- Assemble ---
    intub_extub_df = (
        _first_intubation
        .join(_first_extubation_combined, on="hospitalization_id", how="left")
        .join(_no_extub_classified, on="hospitalization_id", how="left")
        .with_columns(
            pl.when(pl.col("extubation_time").is_not_null())
            .then((pl.col("extubation_time") - pl.col("intubation_time")).dt.total_hours())
            .otherwise(pl.lit(None).cast(pl.Float64))
            .alias("imv_duration_hours")
        )
    )

    n_excluded_no_extub = intub_extub_df.filter(pl.col("extubation_time").is_null()).height

    print(f"Intubations: {intub_extub_df.height:,}")
    print(f"Extubations: {intub_extub_df.filter(pl.col('extubation_time').is_not_null()).height:,}")
    print(f"No extubation detected: {n_excluded_no_extub:,}")
    return intub_extub_df, n_excluded_no_extub


@app.cell
def _(DATA_DIR, FILETYPE, OUTPUT_PHI, TIMEZONE, hosp_final, pl, strip_tz):
    from clifpy import Vitals

    TARGET_VITALS = ["height_cm", "weight_kg", "spo2"]

    _hosp_ids = hosp_final["hospitalization_id"].to_list()
    _vitals = Vitals.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE,
        filters={"hospitalization_id": _hosp_ids}, verbose=False,
    )
    vitals_df = pl.from_pandas(strip_tz(_vitals.df)).filter(
        pl.col("vital_category").str.to_lowercase().is_in(TARGET_VITALS)
    )
    print(f"Vitals loaded (height, weight, spo2): {len(vitals_df):,}")
    for cat in TARGET_VITALS:
        n = vitals_df.filter(pl.col("vital_category").str.to_lowercase() == cat).height
        print(f"  {cat}: {n:,}")

    # BMI: latest weight & height per hospitalization
    _weight = (
        vitals_df
        .filter(pl.col("vital_category").str.to_lowercase() == "weight_kg")
        .sort("recorded_dttm")
        .group_by("hospitalization_id")
        .last()
        .select(["hospitalization_id", pl.col("vital_value").alias("weight_kg")])
    )
    _height = (
        vitals_df
        .filter(pl.col("vital_category").str.to_lowercase() == "height_cm")
        .sort("recorded_dttm")
        .group_by("hospitalization_id")
        .last()
        .select(["hospitalization_id", pl.col("vital_value").alias("height_cm")])
    )
    bmi_df = (
        _weight.join(_height, on="hospitalization_id", how="outer_coalesce")
        .with_columns(
            (pl.col("weight_kg") / (pl.col("height_cm") / 100.0) ** 2).alias("bmi")
        )
    )
    print(f"\nBMI computed: {bmi_df.filter(pl.col('bmi').is_not_null()).height:,} / {bmi_df.height:,} hospitalizations")

    _spo2_df = vitals_df.filter(pl.col("vital_category").str.to_lowercase() == "spo2")
    _spo2_df.write_parquet(str(OUTPUT_PHI / "vitals_spo2_cohort.parquet"))
    print(f"Saved spo2 vitals ({len(_spo2_df):,} rows) to {OUTPUT_PHI / 'vitals_spo2_cohort.parquet'}")
    return bmi_df, vitals_df


@app.cell
def _(
    OUTPUT_PHI,
    adt_df,
    bmi_df,
    first_icu,
    hosp_final,
    intub_extub_df,
    patient_df,
    pl,
    rs_df,
):
    # Re-filter rs_df to only the selected first-ICU-with-IMV window
    _rs_first = (
        rs_df
        .join(
            first_icu.select(["hospitalization_id", "icu_in_dttm", "icu_out_dttm"]),
            on=["hospitalization_id", "icu_in_dttm", "icu_out_dttm"],
            how="semi",
        )
    )
    _rs_counts = (
        _rs_first
        .group_by("hospitalization_id")
        .agg(pl.col("recorded_dttm").count().alias("n_resp_records_first_icu"))
    )

    # Post-ICU location: first ADT row after the ICU window ends
    _post_icu_loc = (
        adt_df
        .join(
            first_icu.select(["hospitalization_id", "icu_out_dttm"]),
            on="hospitalization_id",
            how="inner",
        )
        .filter(pl.col("in_dttm") >= pl.col("icu_out_dttm"))
        .sort(["hospitalization_id", "in_dttm"])
        .unique(subset=["hospitalization_id"], keep="first")
        .select([
            "hospitalization_id",
            pl.col("location_category").alias("post_icu_location"),
        ])
    )

    # Hospital info + ICU location type from first ICU ADT row in the first qualifying ICU stay
    _icu_info = (
        adt_df
        .join(
            first_icu.select(["hospitalization_id", "icu_in_dttm", "icu_out_dttm"]),
            on="hospitalization_id",
            how="inner",
        )
        .filter(
            (pl.col("location_category").str.to_lowercase() == "icu")
            & (pl.col("in_dttm") >= pl.col("icu_in_dttm"))
            & (pl.col("in_dttm") <= pl.col("icu_out_dttm"))
        )
        .sort(["hospitalization_id", "in_dttm"])
        .unique(subset=["hospitalization_id"], keep="first")
        .select(["hospitalization_id", "hospital_id", "hospital_type", "location_type"])
    )

    cohort_df = (
        hosp_final
        .select(["hospitalization_id", "patient_id", "admission_dttm", "discharge_dttm", "age_at_admission", "discharge_category"])
        .join(
            first_icu.select([
                "hospitalization_id",
                "icu_rank",
                pl.col("icu_in_dttm").alias("first_icu_start"),
                pl.col("icu_out_dttm").alias("first_icu_end"),
            ]),
            on="hospitalization_id",
            how="left",
        )
        .join(_rs_counts, on="hospitalization_id", how="left")
        .join(_post_icu_loc, on="hospitalization_id", how="left")
        .join(
            patient_df.select(["patient_id", "sex_category", "race_category", "ethnicity_category", "language_category"]),
            on="patient_id",
            how="left",
        )
        .join(intub_extub_df, on="hospitalization_id", how="left")
        .join(bmi_df, on="hospitalization_id", how="left")
        .join(_icu_info, on="hospitalization_id", how="left")
        .with_columns(
            ((pl.col("discharge_dttm") - pl.col("admission_dttm")).dt.total_hours()).alias("hospital_los_hours"),
            ((pl.col("first_icu_end") - pl.col("first_icu_start")).dt.total_hours()).alias("first_icu_los_hours"),
            pl.col("n_resp_records_first_icu").fill_null(0),
        )
        .with_columns(
            pl.coalesce(["post_icu_location", "discharge_category"]).alias("post_icu_location"),
        )
    )

    # Exclude hospitalizations with no extubation detected before ICU end
    cohort_df = cohort_df.filter(pl.col("extubation_time").is_not_null())

    cohort_df.write_parquet(str(OUTPUT_PHI / "cohort.parquet"))
    print(f"Cohort saved to {OUTPUT_PHI / 'cohort.parquet'}")
    print(f"Shape: {cohort_df.shape}, Columns: {cohort_df.columns}")
    return (cohort_df,)


@app.cell
def _(
    OUTPUT_DIR,
    hosp_24h,
    hosp_adults,
    hosp_after_merge,
    hosp_dated,
    hosp_final,
    hosp_icu,
    json,
    n_excluded_age,
    n_excluded_collar,
    n_excluded_dates,
    n_excluded_no_extub,
    n_excluded_no_icu,
    n_excluded_no_imv,
    n_excluded_null_out,
    n_excluded_short_icu,
    n_excluded_trach,
    n_total,
):
    consort = [
        {
            "step": 0,
            "description": "Total hospitalizations in database",
            "remaining": n_total,
            "excluded": 0,
            "reason": None,
        },
        {
            "step": 1,
            "description": "Adults (age >= 18)",
            "remaining": len(hosp_adults),
            "excluded": n_excluded_age,
            "reason": "Age < 18",
        },
        {
            "step": 2,
            "description": "Admissions within 2018-2024",
            "remaining": len(hosp_dated),
            "excluded": n_excluded_dates,
            "reason": "Admission or discharge outside 2018-2024",
        },
        {
            "step": 3,
            "description": "Has at least one ICU stay",
            "remaining": len(hosp_icu),
            "excluded": n_excluded_no_icu,
            "reason": "No ICU stay",
        },
        {
            "step": 4,
            "description": "ICU stays with valid out_dttm",
            "remaining": len(hosp_after_merge),
            "excluded": n_excluded_null_out,
            "reason": "All ICU stays had null out_dttm",
        },
        {
            "step": 5,
            "description": "At least one ICU stay >= 24 hours",
            "remaining": len(hosp_24h),
            "excluded": n_excluded_short_icu,
            "reason": "No ICU stay >= 24 hours",
        },
        {
            "step": 6,
            "description": "No tracheostomy in ICU stays",
            "remaining": len(hosp_24h) - n_excluded_trach,
            "excluded": n_excluded_trach,
            "reason": "Tracheostomy recorded during ICU stay",
        },
        {
            "step": 7,
            "description": "No trach collar in ICU stays",
            "remaining": len(hosp_24h) - n_excluded_trach - n_excluded_collar,
            "excluded": n_excluded_collar,
            "reason": "Trach collar device during ICU stay",
        },
        {
            "step": 8,
            "description": "Has IMV in an ICU stay >= 24h",
            "remaining": len(hosp_final),
            "excluded": n_excluded_no_imv,
            "reason": "No IMV device in any ICU stay >= 24h",
        },
        {
            "step": 9,
            "description": "Extubation detected before ICU end",
            "remaining": len(hosp_final) - n_excluded_no_extub,
            "excluded": n_excluded_no_extub,
            "reason": "No extubation detected before ICU end",
        },
    ]

    _path = OUTPUT_DIR / "consort.json"
    with open(_path, "w") as _f:
        json.dump(consort, _f, indent=2)
    print(f"CONSORT saved to {_path}")

    from helper import plot_consort

    plot_consort(consort, OUTPUT_DIR / "consort.png")
    print(f"CONSORT diagram saved to {OUTPUT_DIR / 'consort.png'}")
    return (consort,)


@app.cell
def _(consort, pl):
    pl.DataFrame(consort)
    return


@app.cell
def _(
    DATA_DIR,
    FILETYPE,
    OUTPUT_PHI,
    TIMEZONE,
    cohort_df,
    pl,
    strip_tz,
    vitals_df,
):
    from clifpy import MedicationAdminContinuous
    from clifpy.utils.unit_converter import convert_dose_units_by_med_category

    TARGET_MEDS = [
        "norepinephrine", "epinephrine", "phenylephrine", "angiotensin",
        "vasopressin", "dopamine", "dobutamine", "milrinone", "isoproterenol",
        "cisatracurium", "vecuronium", "rocuronium",
        "fentanyl", "propofol", "lorazepam", "midazolam",
        "hydromorphone", "morphine",
    ]

    PREFERRED_UNITS = {
        "norepinephrine": "mcg/kg/min",
        "epinephrine": "mcg/kg/min",
        "phenylephrine": "mcg/kg/min",
        "angiotensin": "ng/kg/min",
        "vasopressin": "u/min",
        "dopamine": "mcg/kg/min",
        "dobutamine": "mcg/kg/min",
        "milrinone": "mcg/kg/min",
        "isoproterenol": "mcg/kg/min",
    }

    # 1) Load continuous medication data
    _mac = MedicationAdminContinuous.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE, verbose=False,
    )
    _mac_pl = pl.from_pandas(strip_tz(_mac.df))
    print(f"Raw medication_admin_continuous records: {len(_mac_pl):,}")

    # 2) Filter to cohort hospitalizations
    _mac_cohort = _mac_pl.join(
        cohort_df.select("hospitalization_id"),
        on="hospitalization_id",
        how="inner",
    )
    print(f"After cohort filter: {len(_mac_cohort):,}")

    # 3) Filter to target med_categories
    _mac_filtered = _mac_cohort.filter(
        pl.col("med_category").str.to_lowercase().is_in(TARGET_MEDS)
    )
    print(f"After med_category filter ({len(TARGET_MEDS)} meds): {len(_mac_filtered):,}")
    print(f"  Med categories found: {sorted(_mac_filtered['med_category'].unique().to_list())}")

    # 4) Vitals for weight-based unit conversion (loaded upstream)
    _vitals_pd = strip_tz(vitals_df.to_pandas())
    print(f"Vitals records for unit conversion: {len(_vitals_pd):,}")

    # 5) Unit conversion (converter works on pandas)
    _mac_pd = _mac_filtered.to_pandas()
    _converted_pd, _counts_pd = convert_dose_units_by_med_category(
        med_df=_mac_pd,
        vitals_df=_vitals_pd,
        preferred_units=PREFERRED_UNITS,
        override=True,
    )

    # 6) Conversion summary
    print(f"\nUnit conversion summary:")
    _status_counts = _converted_pd["_convert_status"].value_counts()
    for status, count in _status_counts.items():
        print(f"  {status}: {count:,}")

    # 7) Convert back to polars
    meds_df = pl.from_pandas(_converted_pd)
    print(f"\nFinal meds_df: {meds_df.shape[0]:,} rows, {meds_df.shape[1]} cols")
    meds_df.write_parquet(str(OUTPUT_PHI / "meds_cohort.parquet"))
    print(f"Saved meds to {OUTPUT_PHI / 'meds_cohort.parquet'}")
    return


@app.cell
def _(
    DATA_DIR,
    FILETYPE,
    OUTPUT_DIR,
    OUTPUT_PHI,
    TIMEZONE,
    cohort_df,
    pl,
    strip_tz,
):
    from clifpy import PatientAssessments

    _ASSESSMENT_CATS = [
        "sat_screen_pass_fail", "sat_screen_performed",
        "sat_delivery_pass_fail", "sat_delivery_performed",
        "sbt_screen_pass_fail", "sbt_screen_performed",
        "sbt_delivery_pass_fail", "sbt_delivery_performed",
        "rass", "gcs_total",
    ]

    _cohort_ids = cohort_df["hospitalization_id"].to_list()
    _pa = PatientAssessments.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE,
        filters={"hospitalization_id": _cohort_ids}, verbose=False,
    )
    assessments_df = pl.from_pandas(strip_tz(_pa.df)).filter(
        pl.col("assessment_category").str.to_lowercase().is_in(_ASSESSMENT_CATS)
    )
    print("Unique values per category (before conversion):")
    _quality_rows = []
    for _cat in _ASSESSMENT_CATS:
        _sub = assessments_df.filter(pl.col("assessment_category").str.to_lowercase() == _cat)
        _num = _sub["numerical_value"].drop_nulls().unique().sort().to_list()
        _catv = _sub["categorical_value"].drop_nulls().unique().sort().to_list()
        print(f"  {_cat} (numerical): {_num}")
        print(f"  {_cat} (categorical): {_catv}")
        _quality_rows.append({"assessment_category": _cat, "numerical_values": ", ".join(str(v) for v in _num), "categorical_values": ", ".join(str(v) for v in _catv)})
    pl.DataFrame(_quality_rows).write_csv(str(OUTPUT_DIR / "assessment_quality.csv"))
    print(f"Saved assessment quality to {OUTPUT_DIR / 'assessment_quality.csv'}")
    # Map categorical pass/fail/yes/no/true/false → 1/0, fill into numerical_value
    _cat_map = {"pass": 1, "yes": 1, "true": 1, "fail": 0, "no": 0, "false": 0}
    _mapped = (
        pl.col("categorical_value")
        .str.to_lowercase()
        .replace(_cat_map)
        .cast(pl.Float64, strict=False)
    )
    assessments_df = assessments_df.with_columns(
        pl.col("numerical_value").fill_null(_mapped).alias("numerical_value"),
    )
    print(f"Patient assessments loaded: {len(assessments_df):,}")
    for _cat in _ASSESSMENT_CATS:
        _n = assessments_df.filter(
            pl.col("assessment_category").str.to_lowercase() == _cat
        ).height
        print(f"  {_cat}: {_n:,}")
    assessments_df.write_parquet(str(OUTPUT_PHI / "assessments_cohort.parquet"))
    print(f"Saved assessments to {OUTPUT_PHI / 'assessments_cohort.parquet'}")
    return


@app.cell
def _(cohort_df):
    cohort_df
    return


if __name__ == "__main__":
    app.run()
