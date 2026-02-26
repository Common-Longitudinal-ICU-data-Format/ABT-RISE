import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import polars as pl
    import json
    from datetime import timedelta
    from pathlib import Path

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
        pl,
        timedelta,
    )


@app.cell
def _(DATA_DIR, FILETYPE, TIMEZONE, pl):
    from clifpy import Hospitalization

    _hosp = Hospitalization.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE, verbose=False,
    )
    hosp_df = pl.from_pandas(_hosp.df)
    n_total = len(hosp_df)
    print(f"Total hospitalizations: {n_total:,}")
    return hosp_df, n_total


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
def _(DATA_DIR, FILETYPE, TIMEZONE, hosp_dated, pl):
    from clifpy import Adt

    _adt = Adt.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE, verbose=False,
    )
    adt_df = pl.from_pandas(_adt.df)

    _icu_hosp_ids = (
        adt_df
        .filter(pl.col("location_category") == "icu")
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
            .is_in(["icu", "procedural"])
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
    icu_stays = _icu_blocks.filter(pl.col("icu_out_dttm").is_not_null())
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
    # Get earliest ICU stay per hospitalization (sort then deduplicate)
    first_icu = (
        icu_stays
        .sort(["hospitalization_id", "icu_in_dttm"])
        .unique(subset=["hospitalization_id"], keep="first")
        .with_columns(
            (pl.col("icu_out_dttm") - pl.col("icu_in_dttm")).alias("icu_los")
        )
    )
    _qualifying = first_icu.filter(pl.col("icu_los") >= timedelta(hours=24))
    hosp_24h = hosp_after_merge.join(
        _qualifying.select("hospitalization_id"),
        on="hospitalization_id",
        how="inner",
    )
    first_icu = _qualifying
    n_excluded_short_icu = len(hosp_after_merge) - len(hosp_24h)
    print(f"Excluded (first ICU < 24h): {n_excluded_short_icu:,}")
    print(f"Remaining: {len(hosp_24h):,}")
    return first_icu, hosp_24h, n_excluded_short_icu


@app.cell
def _(DATA_DIR, FILETYPE, TIMEZONE, first_icu, pl):
    from clifpy import RespiratorySupport

    _rs = RespiratorySupport.from_file(
        data_directory=DATA_DIR, filetype=FILETYPE, timezone=TIMEZONE, verbose=False,
    )
    _rs_all = pl.from_pandas(_rs.df)

    rs_df = (
        _rs_all
        .join(
            first_icu.select(["hospitalization_id", "icu_in_dttm", "icu_out_dttm"]),
            on="hospitalization_id",
            how="inner",
        )
        .filter(
            (pl.col("recorded_dttm") >= pl.col("icu_in_dttm"))
            & (pl.col("recorded_dttm") <= pl.col("icu_out_dttm"))
        )
    )
    print(f"Respiratory support records in first ICU window: {len(rs_df):,}")
    return (rs_df,)


@app.cell
def _(hosp_24h, pl, rs_df):
    _imv_hosp = (
        rs_df
        .filter(pl.col("device_category") == "IMV")
        .select("hospitalization_id")
        .unique()
    )
    hosp_final = hosp_24h.join(_imv_hosp, on="hospitalization_id", how="inner")
    n_excluded_no_imv = len(hosp_24h) - len(hosp_final)
    print(f"Excluded (no IMV in first ICU): {n_excluded_no_imv:,}")
    print(f"Final cohort: {len(hosp_final):,}")
    return hosp_final, n_excluded_no_imv


@app.cell
def _(OUTPUT_PHI, first_icu, hosp_final, pl, rs_df):
    _rs_counts = (
        rs_df
        .group_by("hospitalization_id")
        .agg(pl.col("recorded_dttm").count().alias("n_resp_records_first_icu"))
    )

    cohort_df = (
        hosp_final
        .select(["hospitalization_id", "admission_dttm", "discharge_dttm"])
        .join(
            first_icu.select([
                "hospitalization_id",
                pl.col("icu_in_dttm").alias("first_icu_start"),
                pl.col("icu_out_dttm").alias("first_icu_end"),
            ]),
            on="hospitalization_id",
            how="left",
        )
        .join(_rs_counts, on="hospitalization_id", how="left")
        .with_columns(
            ((pl.col("discharge_dttm") - pl.col("admission_dttm")).dt.total_hours()).alias("hospital_los_hours"),
            ((pl.col("first_icu_end") - pl.col("first_icu_start")).dt.total_hours()).alias("first_icu_los_hours"),
            pl.col("n_resp_records_first_icu").fill_null(0),
        )
    )

    cohort_df.write_parquet(str(OUTPUT_PHI / "cohort.parquet"))
    print(f"Cohort saved to {OUTPUT_PHI / 'cohort.parquet'}")
    print(f"Shape: {cohort_df.shape}, Columns: {cohort_df.columns}")
    return


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
    n_excluded_dates,
    n_excluded_no_icu,
    n_excluded_no_imv,
    n_excluded_null_out,
    n_excluded_short_icu,
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
            "description": "First ICU stay >= 24 hours",
            "remaining": len(hosp_24h),
            "excluded": n_excluded_short_icu,
            "reason": "First ICU stay < 24 hours",
        },
        {
            "step": 6,
            "description": "Has IMV in first ICU stay",
            "remaining": len(hosp_final),
            "excluded": n_excluded_no_imv,
            "reason": "No IMV device in first ICU stay",
        },
    ]

    _path = OUTPUT_DIR / "consort.json"
    with open(_path, "w") as _f:
        json.dump(consort, _f, indent=2)
    print(f"CONSORT saved to {_path}")
    return (consort,)


@app.cell
def _(consort, pl):
    pl.DataFrame(consort)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
