import marimo

__generated_with = "0.20.2"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    import json
    from pathlib import Path

    import polars as pl

    from helper import get_daily_sofa

    # Read site-specific config (site name, paths, etc.)
    _config_path = Path(__file__).parent.parent / "clif_config.json"
    with open(_config_path) as _f:
        _cfg = json.load(_f)

    SITE_NAME = _cfg["site_name"]
    DATA_DIR = _cfg["data_directory"]
    TIMEZONE = _cfg.get("timezone", None)
    SOFA_BATCH_SIZE = int(_cfg.get("sofa_batch_size", 10_000))
    OUTPUT_PHI = Path(__file__).parent.parent / "output_phi"

    print(f"Site: {SITE_NAME}")
    print(f"Data dir: {DATA_DIR}")
    print(f"SOFA cache: {OUTPUT_PHI / 'daily_sofa.parquet'}")
    return DATA_DIR, OUTPUT_PHI, SOFA_BATCH_SIZE, TIMEZONE, get_daily_sofa, pl


@app.cell
def _(OUTPUT_PHI, pl):
    wide = pl.read_parquet(OUTPUT_PHI / "wide_dataset.parquet")
    print(f"Wide dataset: {wide.height:,} rows x {wide.width} cols")

    # Same preprocessing as 03_sat.py / 04_sbt_both.py so the day keys line up
    # (identical casts + sort → identical hosp_id_icu_day set → cache is reused).
    df = wide.with_columns(
        pl.col("recorded_dttm").cast(pl.Datetime("us")),
        pl.col("first_icu_start").cast(pl.Datetime("us")),
        pl.col("first_icu_end").cast(pl.Datetime("us")),
    ).sort("hospitalization_id", "recorded_dttm")
    print(
        f"Preprocessed: {df.height:,} rows, "
        f"{df['hosp_id_icu_day'].n_unique():,} hosp-day keys"
    )
    return (df,)


@app.cell
def _(DATA_DIR, OUTPUT_PHI, SOFA_BATCH_SIZE, TIMEZONE, df, get_daily_sofa, pl):
    # One row per ICU calendar day, worst SOFA computed over [start, end).
    # ── IMPORTANT ──
    # This _day_cohort MUST stay identical to the _day_cohort built in
    # 03_sat.py and 04_sbt_both.py. get_daily_sofa validates its cache by an
    # EXACT hosp_id_icu_day match; if these drift, 03/04 silently recompute
    # SOFA instead of reusing this cache — defeating the whole purpose of this
    # step. Keep the three blocks in sync.
    day_cohort = (
        df.select("hosp_id_icu_day", "hospitalization_id", "icu_day_date")
        .unique(subset=["hosp_id_icu_day"])
        .with_columns(
            pl.col("icu_day_date").alias("start_dttm"),
            (pl.col("icu_day_date") + pl.duration(days=1)).alias("end_dttm"),
        )
    )
    print(f"Day cohort: {day_cohort.height:,} unique ICU days")

    # refresh=True: this step's job is to (re)build the shared cache. Run it
    # after 02_wide_dataset.py and before 03_sat.py / 04_sbt_both.py.
    daily_sofa = get_daily_sofa(
        data_directory=DATA_DIR,
        day_cohort=day_cohort,
        cache_path=OUTPUT_PHI / "daily_sofa.parquet",
        source_path=OUTPUT_PHI / "wide_dataset.parquet",
        timezone=TIMEZONE,
        batch_size=SOFA_BATCH_SIZE,
        refresh=True,
    )

    _total = daily_sofa["sofa_total"]
    _nn = _total.drop_nulls()
    print(f"\nDaily SOFA cached: {daily_sofa.height:,} days")
    print(
        f"  with SOFA: {_nn.len():,}  |  null: {daily_sofa.height - _nn.len():,}"
    )
    if _nn.len():
        print(
            f"  sofa_total: min={_nn.min()}, median={_nn.median()}, max={_nn.max()}"
        )
    return daily_sofa, day_cohort


if __name__ == "__main__":
    app.run()
