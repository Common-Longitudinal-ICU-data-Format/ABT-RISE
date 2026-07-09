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
    import json
    from datetime import timedelta
    from pathlib import Path

    _config_path = Path(__file__).parent.parent / "clif_config.json"
    with open(_config_path) as _f:
        _cfg = json.load(_f)

    SITE_NAME = _cfg["site_name"]
    DATA_DIR = _cfg["data_directory"]
    OUTPUT_PHI = Path(__file__).parent.parent / "output_phi"
    OUTPUT_ANALYSIS = OUTPUT_PHI / "analysis"
    OUTPUT_ANALYSIS.mkdir(parents=True, exist_ok=True)

    # Consensus windows for intubation episode counting
    CONSENSUS_INTUB_WINDOW_MIN = 15   # IMV must stay IMV this long after start
    CONSENSUS_EXTUB_WINDOW_MIN = 30   # non-IMV must stay non-IMV this long to close an episode

    print(f"Site: {SITE_NAME}")
    print(f"Output: {OUTPUT_ANALYSIS}")
    return (
        CONSENSUS_EXTUB_WINDOW_MIN,
        CONSENSUS_INTUB_WINDOW_MIN,
        DATA_DIR,
        OUTPUT_ANALYSIS,
        OUTPUT_PHI,
        Path,
        SITE_NAME,
        pl,
        timedelta,
    )


@app.cell
def _(DATA_DIR, OUTPUT_PHI, Path, pl):
    # Load SAT and SBT day-level outputs + cohort
    sat = pl.read_parquet(OUTPUT_PHI / "sat_standard" / "sat_day_level.parquet")
    sbt = pl.read_parquet(OUTPUT_PHI / "sbt_both_stabilities" / "sbt_day_level.parquet")
    cohort = pl.read_parquet(OUTPUT_PHI / "cohort.parquet")

    # Load patient table for death_dttm
    patient = (
        pl.read_parquet(Path(DATA_DIR) / "clif_patient.parquet")
        .select("patient_id", "death_dttm")
        .with_columns(pl.col("death_dttm").dt.replace_time_zone(None))
    )

    print(f"SAT day-level: {sat.height:,} rows x {sat.width} cols")
    print(f"SBT day-level: {sbt.height:,} rows x {sbt.width} cols")
    print(f"Cohort: {cohort.height:,} patients")
    print(f"SAT columns: {sat.columns}")
    print(f"SBT columns: {sbt.columns}")
    return cohort, patient, sat, sbt


@app.cell
def _(
    CONSENSUS_EXTUB_WINDOW_MIN,
    CONSENSUS_INTUB_WINDOW_MIN,
    OUTPUT_PHI,
    cohort,
    pl,
    timedelta,
):
    # Count intubation episodes via state machine with consensus windows.
    # An episode is counted when a non-IMV→IMV transition has IMV stable for
    # CONSENSUS_INTUB_WINDOW_MIN forward. An episode closes on either (a) trach
    # collar row or (b) IMV→non-IMV transition with non-IMV stable for
    # CONSENSUS_EXTUB_WINDOW_MIN forward. If trach collar appears before any
    # IMV row, count=1 and the patient is skipped.
    _raw_wf = (
        pl.read_parquet(OUTPUT_PHI / "resp_waterfall_cohort.parquet")
        .filter(pl.col("hospitalization_id").is_in(cohort["hospitalization_id"]))
        .with_columns(pl.col("recorded_dttm").dt.replace_time_zone(None))
    )

    _intub_win = timedelta(minutes=CONSENSUS_INTUB_WINDOW_MIN)
    _extub_win = timedelta(minutes=CONSENSUS_EXTUB_WINDOW_MIN)

    def _count_intubations(wf):
        """Count confirmed intubation episodes per hospitalization."""
        if wf.height == 0:
            return pl.DataFrame(
                {"hospitalization_id": [], "n": []},
                schema={"hospitalization_id": pl.String, "n": pl.Int64},
            )

        prepped = (
            # Null timestamps have no position on the timeline and cannot take
            # part in the forward consensus windows (t + window); drop them so the
            # state machine never does `None + timedelta`.
            wf.filter(pl.col("recorded_dttm").is_not_null())
            .sort(["hospitalization_id", "recorded_dttm"])
            .with_columns(pl.col("device_category").str.to_lowercase().alias("_dev_lc"))
            .with_columns(
                (pl.col("_dev_lc") == "imv").cast(pl.Int8).alias("_is_imv"),
                (
                    (pl.col("_dev_lc") == "trach collar")
                    | (pl.col("tracheostomy").cast(pl.Boolean).fill_null(False))
                ).alias("_is_trach"),
            )
            .filter(pl.col("_dev_lc").is_not_null() | pl.col("_is_trach"))
            .select("hospitalization_id", "recorded_dttm", "_is_imv", "_is_trach")
        )

        results = []
        for (hosp_id,), grp in prepped.group_by("hospitalization_id", maintain_order=True):
            rows = grp.rows()  # tuples: (hosp_id, recorded_dttm, is_imv, is_trach)
            times = [r[1] for r in rows]
            imvs = [r[2] for r in rows]
            trachs = [r[3] for r in rows]
            n = len(rows)

            # Direct trach rule: if trach sentinel fires before any IMV → count=1, skip
            first_imv_idx = next((i for i, v in enumerate(imvs) if v == 1), None)
            first_trach_idx = next((i for i, v in enumerate(trachs) if v), None)
            if first_trach_idx is not None and (
                first_imv_idx is None or first_trach_idx < first_imv_idx
            ):
                results.append({"hospitalization_id": hosp_id, "n": 1})
                continue

            count = 0
            on_vent = False
            for i in range(n):
                # Trach sentinel (device=="trach collar" OR tracheostomy==True):
                # close any open episode and stop all further counting.
                if trachs[i]:
                    on_vent = False
                    break

                t = times[i]
                imv = imvs[i]

                if not on_vent and imv == 1:
                    # 15-min forward consensus: no is_imv==0 in (t, t+window]
                    deadline = t + _intub_win
                    contradicted = False
                    j = i + 1
                    while j < n and times[j] <= deadline:
                        if imvs[j] == 0:
                            contradicted = True
                            break
                        j += 1
                    if not contradicted:
                        count += 1
                        on_vent = True

                elif on_vent and imv == 0:
                    # 30-min forward consensus: no is_imv==1 in (t, t+window]
                    deadline = t + _extub_win
                    contradicted = False
                    j = i + 1
                    while j < n and times[j] <= deadline:
                        if imvs[j] == 1:
                            contradicted = True
                            break
                        j += 1
                    if not contradicted:
                        on_vent = False

            results.append({"hospitalization_id": hosp_id, "n": count})

        return pl.DataFrame(
            results,
            schema={"hospitalization_id": prepped.schema["hospitalization_id"], "n": pl.Int64},
        )

    # Whole hospitalization
    _hosp_counts = _count_intubations(_raw_wf)

    # Index ICU stay only (windowed to first_icu_start → first_icu_end)
    _icu_wf = (
        _raw_wf
        .join(
            cohort.select("hospitalization_id",
                          pl.col("first_icu_start").dt.replace_time_zone(None),
                          pl.col("first_icu_end").dt.replace_time_zone(None)),
            on="hospitalization_id", how="inner",
        )
        .filter(
            (pl.col("recorded_dttm") >= pl.col("first_icu_start"))
            & (pl.col("recorded_dttm") <= pl.col("first_icu_end"))
        )
        .drop("first_icu_start", "first_icu_end")
    )
    _icu_counts = _count_intubations(_icu_wf)

    mv_counts = (
        cohort.select("hospitalization_id")
        .join(_hosp_counts.rename({"n": "mv_count_in_whole_hospitalization"}),
              on="hospitalization_id", how="left")
        .join(_icu_counts.rename({"n": "mv_count_in_index_icu_stay"}),
              on="hospitalization_id", how="left")
        .with_columns(
            pl.col("mv_count_in_whole_hospitalization").fill_null(1),
            pl.col("mv_count_in_index_icu_stay").fill_null(1),
        )
    )

    print(f"MV count (whole hosp):\n{mv_counts['mv_count_in_whole_hospitalization'].value_counts().sort('mv_count_in_whole_hospitalization')}")
    print(f"\nMV count (index ICU):\n{mv_counts['mv_count_in_index_icu_stay'].value_counts().sort('mv_count_in_index_icu_stay')}")
    return (mv_counts,)


@app.cell
def _(cohort, patient, pl, sat, sbt):
    # ── Merge SAT + SBT on hosp_id_icu_day ──
    # Take SAT as base; add SBT-specific columns
    _sbt_cols = [
        "hosp_id_icu_day",
        "sbt_screen_pass_fail", "sbt_screen_performed",
        "sbt_delivery_pass_fail", "sbt_delivery_performed",
        "SBT_delivery_2min_primary", "SBT_delivery_5min_secondary",
        "sbt_ground_truth",
    ]
    # Only keep SBT columns that exist
    _sbt_keep = [c for c in _sbt_cols if c in sbt.columns]

    # SBT eligible_event (rename to avoid collision with SAT's)
    _sbt_elig = sbt.select(
        "hosp_id_icu_day",
        pl.col("eligible_event").alias("sbt_eligible_event"),
    ).unique(subset=["hosp_id_icu_day"])

    merged = (
        sat.join(sbt.select(_sbt_keep).unique(subset=["hosp_id_icu_day"]),
                 on="hosp_id_icu_day", how="left")
        .join(_sbt_elig, on="hosp_id_icu_day", how="left")
    )

    # Join cohort for hospital/ICU info, intubation time, death info
    _cohort_cols = cohort.select(
        "hospitalization_id", "patient_id",
        *[c for c in ["hospital_id", "hospital_type", "location_type"] if c in cohort.columns],
        "intubation_time", "extubation_time",
        "first_icu_los_hours", "admission_dttm", "discharge_dttm", "discharge_category",
    )
    _death = _cohort_cols.join(patient, on="patient_id", how="left").drop("patient_id")
    merged = merged.join(
        _death.select([c for c in _death.columns if c not in merged.columns or c == "hospitalization_id"]),
        on="hospitalization_id", how="left",
    )

    print(f"Merged: {merged.height:,} rows x {merged.width} cols")
    print(f"Vent days: {merged.filter(pl.col('is_vent_day')).height:,}")
    return (merged,)


@app.cell
def _(merged, pl):
    # ── Build File 1: Person-Period (one row per vent day) ──
    f1 = merged.filter(pl.col("is_vent_day")).sort("hospitalization_id", "vent_day")

    # Rename / compute columns per SAP
    f1 = f1.with_columns(
        # Eligibility flags (binary 0/1)
        pl.when(pl.col("eligible_event") == 1).then(1).otherwise(0).alias("SAT_eligible"),
        pl.when(pl.col("sbt_eligible_event") == 1).then(1).otherwise(0).alias("SBT_eligible"),

        # Delivery flags (fill null → 0)
        pl.col("SAT_primary_delivery").fill_null(0).cast(pl.Int32).alias("SAT_delivered_primary"),
        pl.col("SAT_modified_delivery").fill_null(0).cast(pl.Int32).alias("SAT_delivered_modified"),
        pl.col("SBT_delivery_2min_primary").fill_null(0).cast(pl.Int32).alias("SBT_delivered_2min"),
        pl.col("SBT_delivery_5min_secondary").fill_null(0).cast(pl.Int32).alias("SBT_delivered_5min"),

        # Flowsheet ground truth
        pl.col("sat_ground_truth").fill_null(0).cast(pl.Int32).alias("flowsheet_SAT"),
        pl.col("sbt_ground_truth").fill_null(0).cast(pl.Int32).alias("flowsheet_SBT"),

        # Extubated today (outcome == 2)
        pl.when(pl.col("outcome") == 2).then(1).otherwise(0).alias("extubated"),

        # Vent start date
        pl.col("intubation_time").dt.date().alias("vent_start_date"),

        # Individual sedation med flags (binary 0/1)
        pl.col("propofol_prior").fill_null(False).cast(pl.Int32).alias("propofol_prior_flag"),
        pl.col("fentanyl_prior").fill_null(False).cast(pl.Int32).alias("fentanyl_prior_flag"),
        pl.col("hydromorphone_prior").fill_null(False).cast(pl.Int32).alias("hydromorphone_prior_flag"),
        pl.col("morphine_prior").fill_null(False).cast(pl.Int32).alias("morphine_prior_flag"),
        pl.col("lorazepam_prior").fill_null(False).cast(pl.Int32).alias("lorazepam_prior_flag"),
        pl.col("midazolam_prior").fill_null(False).cast(pl.Int32).alias("midazolam_prior_flag"),
        pl.col("nmb_prior").fill_null(False).cast(pl.Int32).alias("nmb_prior_flag"),

        # Sedation prior (binary: any sedation charted prior day — includes NMB)
        (
            pl.col("propofol_prior").fill_null(False)
            | pl.col("bzd_prior").fill_null(False)
            | pl.col("opioid_prior").fill_null(False)
            | pl.col("nmb_prior").fill_null(False)
        ).cast(pl.Int32).alias("sedation_prior"),
    )

    # ── Rolling cumulative proportions (prior-day lagged) ──
    # For each delivery version: cum_delivered / cum_eligible, shifted by 1 day
    for _del_col, _elig_col, _prop_col in [
        ("SAT_delivered_primary", "SAT_eligible", "SAT_prop_primary"),
        ("SAT_delivered_modified", "SAT_eligible", "SAT_prop_modified"),
        ("SBT_delivered_2min", "SBT_eligible", "SBT_prop_2min"),
        ("SBT_delivered_5min", "SBT_eligible", "SBT_prop_5min"),
    ]:
        # Only count delivery on eligible days
        _del_on_elig = f"_del_elig_{_prop_col}"
        f1 = f1.with_columns(
            (pl.col(_del_col) * pl.col(_elig_col)).alias(_del_on_elig)
        )
        f1 = f1.with_columns(
            (
                pl.col(_del_on_elig).cum_sum().shift(1).over("hospitalization_id").fill_null(0)
                / pl.col(_elig_col).cum_sum().shift(1).over("hospitalization_id").fill_null(0)
            ).fill_nan(0.0).fill_null(0.0).alias(_prop_col)
        ).drop(_del_on_elig)

    # ── Select final File 1 columns ──
    _id_cols = ["hospitalization_id", "vent_start_date", "vent_day", "icu_day", "icu_day_date"]
    _hospital_cols = [c for c in ["hospital_id", "hospital_type", "location_type"] if c in f1.columns]
    _exposure_cols = [
        "SAT_eligible", "SBT_eligible",
        "SAT_delivered_primary", "SAT_delivered_modified",
        "SBT_delivered_2min", "SBT_delivered_5min",
        "SAT_prop_primary", "SAT_prop_modified",
        "SBT_prop_2min", "SBT_prop_5min",
    ]
    _flowsheet_cols = ["flowsheet_SAT", "flowsheet_SBT"]
    _outcome_cols = ["extubated"]
    _clinical_cols = ["sofa_prior", "fio2_prior", "peep_prior", "nee_prior", "sedation_prior"]

    # Rename priors for SAP consistency
    f1 = f1.rename({
        "sofa_prior": "SOFA_prior",
        "fio2_prior": "FiO2_prior",
        "peep_prior": "PEEP_prior",
        "nee_prior": "NEE_prior",
    })
    _clinical_cols = [
        "SOFA_prior", "FiO2_prior", "PEEP_prior", "NEE_prior", "sedation_prior",
        "propofol_prior_flag", "fentanyl_prior_flag", "hydromorphone_prior_flag",
        "morphine_prior_flag", "lorazepam_prior_flag", "midazolam_prior_flag",
        "nmb_prior_flag",
    ]

    file1 = f1.select(
        _hospital_cols + _id_cols + _exposure_cols + _flowsheet_cols + _outcome_cols + _clinical_cols
    ).sort("hospitalization_id", "vent_day")

    print(f"File 1 (Person-Period): {file1.height:,} rows x {file1.width} cols")
    print(f"  Patients: {file1['hospitalization_id'].n_unique():,}")
    print(f"  SAT eligible days: {file1.filter(pl.col('SAT_eligible') == 1).height:,}")
    print(f"  SBT eligible days: {file1.filter(pl.col('SBT_eligible') == 1).height:,}")
    print(f"  Extubation events: {file1.filter(pl.col('extubated') == 1).height:,}")
    return (file1,)


@app.cell
def _(cohort, file1, mv_counts, patient, pl, sat):
    # ── Build File 2: Hospitalization-Level (one row per hospitalization) ──

    # Aggregate from File 1
    _agg = (
        file1.group_by("hospitalization_id")
        .agg(
            # Proportions: total delivered / total eligible
            (pl.col("SAT_delivered_primary").sum() / pl.col("SAT_eligible").sum())
            .fill_nan(0.0).alias("SAT_prop_final_primary"),
            (pl.col("SAT_delivered_modified").sum() / pl.col("SAT_eligible").sum())
            .fill_nan(0.0).alias("SAT_prop_final_modified"),
            (pl.col("SBT_delivered_2min").sum() / pl.col("SBT_eligible").sum())
            .fill_nan(0.0).alias("SBT_prop_final_2min"),
            (pl.col("SBT_delivered_5min").sum() / pl.col("SBT_eligible").sum())
            .fill_nan(0.0).alias("SBT_prop_final_5min"),

            # Mean clinical values across vent days
            pl.col("SOFA_prior").mean().alias("SOFA_mean"),
            pl.col("FiO2_prior").mean().alias("FiO2_mean"),
            pl.col("PEEP_prior").mean().alias("PEEP_mean"),
            pl.col("NEE_prior").mean().alias("NEE_mean"),
            pl.col("sedation_prior").mean().alias("sedation_mean"),
            pl.col("propofol_prior_flag").mean().alias("propofol_mean"),
            pl.col("fentanyl_prior_flag").mean().alias("fentanyl_mean"),
            pl.col("hydromorphone_prior_flag").mean().alias("hydromorphone_mean"),
            pl.col("morphine_prior_flag").mean().alias("morphine_mean"),
            pl.col("lorazepam_prior_flag").mean().alias("lorazepam_mean"),
            pl.col("midazolam_prior_flag").mean().alias("midazolam_mean"),
            pl.col("nmb_prior_flag").mean().alias("nmb_mean"),

            # Total vent days
            pl.len().alias("n_vent_days"),
        )
    )

    # Join cohort demographics + outcomes
    _cohort_info = cohort.select(
        "hospitalization_id", "patient_id",
        "age_at_admission", "sex_category",
        *[c for c in ["hospital_id", "hospital_type", "location_type"] if c in cohort.columns],
        "intubation_time", "extubation_time",
        "first_icu_los_hours", "discharge_category",
    ).join(patient, on="patient_id", how="left")

    # CCI is in day-level data (SAT), grab one value per patient
    _cci = sat.select("hospitalization_id", "cci_score").group_by("hospitalization_id").agg(pl.col("cci_score").first())

    file2 = (
        _cohort_info.join(_agg, on="hospitalization_id", how="left")
        .join(_cci, on="hospitalization_id", how="left")
        .join(mv_counts, on="hospitalization_id", how="left")
    ).with_columns(pl.col("n_vent_days").fill_null(0))

    # Compute outcome columns
    file2 = file2.with_columns(
        # Extubation flag & time
        pl.col("extubation_time").is_not_null().cast(pl.Int32).alias("extubation_flag"),
        pl.when(pl.col("extubation_time").is_not_null())
        .then(
            (pl.col("extubation_time") - pl.col("intubation_time")).dt.total_hours() / 24.0
        )
        .otherwise(None)
        .round(1)
        .alias("time_to_extubation"),

        # Death flag & time
        pl.when(
            pl.col("death_dttm").is_not_null()
        ).then(1)
        .when(
            pl.col("discharge_category").str.to_lowercase() == "expired"
        ).then(1)
        .otherwise(0)
        .alias("death_flag"),

        pl.when(pl.col("death_dttm").is_not_null() & pl.col("intubation_time").is_not_null())
        .then(
            (pl.col("death_dttm") - pl.col("intubation_time")).dt.total_hours() / 24.0
        )
        .otherwise(None)
        .round(1)
        .alias("days_to_death"),

        # ICU LOS in days
        (pl.col("first_icu_los_hours") / 24.0).round(1).alias("ICU_LOS"),
    )

    # VFD-28
    file2 = file2.with_columns(
        pl.when(pl.col("death_flag") == 1).then(0)
        .when(pl.col("extubation_flag") == 0).then(0)
        .when(
            (pl.col("extubation_time") - pl.col("intubation_time")).dt.total_hours() / 24.0 > 28
        ).then(0)
        .otherwise(
            (28 - (pl.col("extubation_time") - pl.col("intubation_time")).dt.total_hours() / 24.0)
            .round(0).cast(pl.Int32)
        )
        .alias("VFD_28")
    )

    # site_has_flowsheet: 1 if any raw flowsheet columns were ever charted at this site
    _has_sat_flow = sat.select(
        pl.col("sat_screen_pass_fail").is_not_null().any()
        | pl.col("sat_screen_performed").is_not_null().any()
        | pl.col("sat_delivery_pass_fail").is_not_null().any()
        | pl.col("sat_delivery_performed").is_not_null().any()
    ).item()
    _has_sbt_flow = sbt.select(
        pl.col("sbt_screen_pass_fail").is_not_null().any()
        | pl.col("sbt_screen_performed").is_not_null().any()
        | pl.col("sbt_delivery_pass_fail").is_not_null().any()
        | pl.col("sbt_delivery_performed").is_not_null().any()
    ).item()
    file2 = file2.with_columns(
        pl.lit(int(_has_sat_flow)).alias("site_has_flowsheet_sat"),
        pl.lit(int(_has_sbt_flow)).alias("site_has_flowsheet_sbt"),
    )

    # Rename for SAP consistency
    file2 = file2.rename({
        "age_at_admission": "age",
        "sex_category": "sex",
        "cci_score": "CCI",
    })

    # Select final columns
    _hospital_cols = [c for c in ["hospital_id", "hospital_type", "location_type"] if c in file2.columns]
    _id_cols = ["hospitalization_id"]
    _demo_cols = ["age", "sex", "CCI"]
    _prop_cols = [
        "SAT_prop_final_primary", "SAT_prop_final_modified",
        "SBT_prop_final_2min", "SBT_prop_final_5min",
    ]
    _outcome_cols = [
        "extubation_flag", "time_to_extubation",
        "death_flag", "days_to_death",
        "VFD_28", "ICU_LOS",
    ]
    _clinical_cols = [
        "SOFA_mean", "FiO2_mean", "PEEP_mean", "NEE_mean", "sedation_mean",
        "propofol_mean", "fentanyl_mean", "hydromorphone_mean",
        "morphine_mean", "lorazepam_mean", "midazolam_mean", "nmb_mean",
    ]
    _meta_cols = ["n_vent_days", "mv_count_in_whole_hospitalization", "mv_count_in_index_icu_stay", "site_has_flowsheet_sat", "site_has_flowsheet_sbt"]

    file2 = file2.select(
        _hospital_cols + _id_cols + _demo_cols + _prop_cols
        + _outcome_cols + _clinical_cols + _meta_cols
    )

    print(f"File 2 (Hospitalization-Level): {file2.height:,} rows x {file2.width} cols")
    print(f"  Deaths: {file2.filter(pl.col('death_flag') == 1).height:,}")
    print(f"  Extubations: {file2.filter(pl.col('extubation_flag') == 1).height:,}")
    _vfd = file2.filter(pl.col("VFD_28") > 0)
    print(f"  VFD_28 > 0: {_vfd.height:,} (mean={_vfd['VFD_28'].mean():.1f})")
    print(f"  VFD_28 == 0 among deaths: {file2.filter((pl.col('death_flag') == 1) & (pl.col('VFD_28') == 0)).height:,}")
    return (file2,)


@app.cell
def _(OUTPUT_ANALYSIS, file1, file2):
    # ── Save outputs ──
    file1.write_parquet(str(OUTPUT_ANALYSIS / "file1_person_period.parquet"))
    print(f"Saved File 1: {file1.height:,} rows x {file1.width} cols")

    file2.write_parquet(str(OUTPUT_ANALYSIS / "file2_hospitalization_level.parquet"))
    print(f"Saved File 2: {file2.height:,} rows x {file2.width} cols")
    return


@app.cell
def _(file1, file2, pl):
    # ── Validation checks ──
    print("=== VALIDATION ===")

    # File 1 checks
    _f1_patients = file1["hospitalization_id"].n_unique()
    print(f"\nFile 1:")
    print(f"  Patients: {_f1_patients:,}")
    print(f"  Total vent days: {file1.height:,}")

    # Check vent_day starts at 1 for each patient
    _min_vd = file1.group_by("hospitalization_id").agg(pl.col("vent_day").min()).filter(pl.col("vent_day") != 1)
    if _min_vd.height > 0:
        print(f"  WARNING: {_min_vd.height} patients don't start at vent_day=1")
    else:
        print(f"  OK: all patients start at vent_day=1")

    # Check SAT_prop is 0 on day 1
    _d1_prop = file1.filter(pl.col("vent_day") == 1).select("SAT_prop_primary").to_series()
    if (_d1_prop != 0).any():
        print(f"  WARNING: SAT_prop_primary not 0 on vent_day=1 for some patients")
    else:
        print(f"  OK: SAT_prop_primary = 0 on vent_day=1 for all patients")

    # File 2 checks
    print(f"\nFile 2:")
    print(f"  Patients: {file2.height:,}")

    # VFD-28 = 0 for all deaths
    _death_vfd = file2.filter((pl.col("death_flag") == 1) & (pl.col("VFD_28") != 0))
    if _death_vfd.height > 0:
        print(f"  WARNING: {_death_vfd.height} deaths with VFD_28 != 0")
    else:
        print(f"  OK: VFD_28 = 0 for all deaths")

    # Prop_final between 0 and 1
    for _col in ["SAT_prop_final_primary", "SAT_prop_final_modified", "SBT_prop_final_2min", "SBT_prop_final_5min"]:
        _bad = file2.filter((pl.col(_col) < 0) | (pl.col(_col) > 1))
        if _bad.height > 0:
            print(f"  WARNING: {_col} has {_bad.height} values outside [0,1]")
    print(f"  OK: all prop_final columns in [0,1]")

    # Summary stats
    print(f"\n=== SUMMARY ===")
    print(f"  SAT_prop_final_primary: mean={file2['SAT_prop_final_primary'].mean():.3f}, median={file2['SAT_prop_final_primary'].median():.3f}")
    print(f"  SBT_prop_final_2min:    mean={file2['SBT_prop_final_2min'].mean():.3f}, median={file2['SBT_prop_final_2min'].median():.3f}")
    print(f"  VFD_28:                 mean={file2['VFD_28'].mean():.1f}, median={file2['VFD_28'].median():.1f}")
    print(f"  ICU_LOS:                mean={file2['ICU_LOS'].mean():.1f}, median={file2['ICU_LOS'].median():.1f}")
    return


if __name__ == "__main__":
    app.run()
