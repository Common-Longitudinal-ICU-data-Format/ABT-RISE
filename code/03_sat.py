import marimo

__generated_with = "0.20.2"
app = marimo.App(width="full")


@app.cell
def _():
    import polars as pl
    import pandas as pd
    import numpy as np
    import json
    from pathlib import Path
    from tqdm import tqdm

    _config_path = Path(__file__).parent.parent / "clif_config.json"
    with open(_config_path) as _f:
        _cfg = json.load(_f)

    SITE_NAME = _cfg["site_name"]
    OUTPUT_PHI = Path(__file__).parent.parent / "output_phi"
    OUTPUT_SAT = OUTPUT_PHI / "sat_standard"
    OUTPUT_SAT.mkdir(parents=True, exist_ok=True)

    print(f"Site: {SITE_NAME}")
    print(f"Output: {OUTPUT_SAT}")
    return OUTPUT_PHI, OUTPUT_SAT, SITE_NAME, np, pd, pl, tqdm


@app.cell
def _(OUTPUT_PHI, pl):
    wide = pl.read_parquet(OUTPUT_PHI / "wide_dataset.parquet")
    print(f"Wide dataset: {wide.height:,} rows x {wide.width} cols")

    cohort = pl.read_parquet(OUTPUT_PHI / "cohort.parquet")
    print(f"Cohort: {cohort.height:,} hospitalizations")
    return cohort, wide


@app.cell
def _(pd, wide):
    # Constants from definitions_source_of_truth.py
    VENT_DAY_ANCHOR_HOUR = 6
    MAX_FFILL_OBSERVATIONS = 6

    df = wide.to_pandas()
    df["recorded_dttm"] = pd.to_datetime(df["recorded_dttm"])
    df["first_icu_start"] = pd.to_datetime(df["first_icu_start"])
    df["first_icu_end"] = pd.to_datetime(df["first_icu_end"])

    # 06:00-anchored vent-day date: shift back 6h then normalize
    df["vent_day_date"] = (df["recorded_dttm"] - pd.Timedelta(hours=VENT_DAY_ANCHOR_HOUR)).dt.normalize()
    df["hosp_id_day_key"] = df["hospitalization_id"].astype(str) + "_" + df["vent_day_date"].dt.strftime("%Y-%m-%d")

    df = df.sort_values(["hospitalization_id", "recorded_dttm"]).reset_index(drop=True)
    print(f"Preprocessed: {len(df):,} rows, {df['hosp_id_day_key'].nunique():,} hosp-day keys")
    return (df,)


@app.cell
def _(pd, tqdm):
    def process_cohort(df: pd.DataFrame) -> pd.DataFrame:
        """Identify encounter-days where SAT eligibility conditions are met
        for at least 4 continuous hours in the 22:00-06:00 window.

        The 22:00-06:00 window spans two 06:00-anchored vent-day groups,
        so we group by hospitalization_id only and iterate over each
        unique vent_day_date, pulling the window data from all rows for
        that hospitalization.
        """
        df = df.sort_values(["hospitalization_id", "recorded_dttm"]).reset_index(drop=True)
        df = df.copy()

        df["all_conditions_check"] = (
            (df["device_category"].str.lower() == "imv")
            & (df["min_sedation_dose_2"] > 0)
            & (df["max_paralytics"] <= 0)
        ).astype(int)

        # Restrict to days with at least one IMV observation
        vented_days = df[df["device_category"].str.lower() == "imv"]["hosp_id_day_key"].unique()
        df = df[df["hosp_id_day_key"].isin(vented_days)]

        result = []
        for hosp_id, hosp_df in tqdm(df.groupby("hospitalization_id"), desc="SAT eligibility check"):
            hosp_df = hosp_df.sort_values("recorded_dttm")
            times = hosp_df["recorded_dttm"].values

            for date in hosp_df["vent_day_date"].unique():
                # Window: previous day 22:00 to current day 06:00
                start_time = date + pd.Timedelta(hours=22) - pd.Timedelta(days=1)
                end_time = date + pd.Timedelta(hours=6)
                mask = (times >= start_time) & (times <= end_time)
                window_df = hosp_df.loc[mask]

                if window_df.empty or not window_df["all_conditions_check"].any():
                    continue

                window_df = window_df.copy()
                window_df["condition_met_group"] = (
                    (window_df["all_conditions_check"] != window_df["all_conditions_check"].shift()).cumsum()
                )
                valid_segments = window_df[window_df["all_conditions_check"] == 1].groupby("condition_met_group")
                for _, segment in valid_segments:
                    segment = segment.sort_values("recorded_dttm").copy()
                    segment["duration"] = segment["recorded_dttm"].diff().fillna(pd.Timedelta(seconds=0))
                    segment["cumulative_duration"] = segment["duration"].cumsum()
                    if segment["cumulative_duration"].iloc[-1] >= pd.Timedelta(hours=4):
                        event_time_at_4_hours = (
                            segment[segment["cumulative_duration"] >= pd.Timedelta(hours=4)]
                            .iloc[0]["recorded_dttm"]
                        )
                        result.append({
                            "hospitalization_id": hosp_id,
                            "current_day_key": date,
                            "event_time_at_4_hours": event_time_at_4_hours,
                        })
                        break

        if not result:
            return pd.DataFrame(columns=["hospitalization_id", "current_day_key", "event_time_at_4_hours"])
        return pd.DataFrame(result)

    return (process_cohort,)


@app.cell
def _(df, np, process_cohort):
    result_df = process_cohort(df)
    print(f"Encounter-days with ≥4h eligibility: {len(result_df):,}")

    cohort_work = df.copy()
    cohort_elig = cohort_work.merge(
        result_df[["hospitalization_id", "current_day_key", "event_time_at_4_hours"]],
        how="left",
        left_on=["hospitalization_id", "vent_day_date"],
        right_on=["hospitalization_id", "current_day_key"],
    )

    _mask_valid = cohort_elig["event_time_at_4_hours"].notna()
    _mask_after = cohort_elig["recorded_dttm"] >= cohort_elig["event_time_at_4_hours"]
    _eligible_rows = cohort_elig[_mask_valid & _mask_after].copy()
    _first_eligible_idx = _eligible_rows.groupby(["hospitalization_id", "vent_day_date"])["recorded_dttm"].idxmin()

    cohort_elig["eligible_event"] = np.nan
    cohort_elig.loc[_first_eligible_idx, "eligible_event"] = 1

    # Remove flag on last event per hospitalization (edge artefact)
    _last_idxs = cohort_elig.groupby("hospitalization_id")["recorded_dttm"].idxmax()
    cohort_elig.loc[_last_idxs, "eligible_event"] = np.nan

    _eligible_days = cohort_elig.loc[cohort_elig["eligible_event"] == 1, "hosp_id_day_key"].unique()
    cohort_elig["on_vent_and_sedation"] = cohort_elig["hosp_id_day_key"].isin(_eligible_days).astype(int)
    cohort_elig = cohort_elig.drop(columns=["current_day_key", "event_time_at_4_hours"], errors="ignore")

    print(f"Eligible hosp-day keys: {cohort_elig[cohort_elig['on_vent_and_sedation'] == 1]['hosp_id_day_key'].nunique():,}")
    return (cohort_elig,)


@app.cell
def _(cohort_elig, np, pd, tqdm):
    _MED_COLS = ["fentanyl", "propofol", "lorazepam", "midazolam", "hydromorphone", "morphine"]
    _MED_COLS2 = ["propofol", "lorazepam", "midazolam"]
    _FLAG_COLS = [
        "SAT_EHR_delivery", "SAT_modified_delivery", "SAT_rass_nonneg_30",
        "SAT_med_halved_rass_pos", "SAT_no_meds_rass_pos_45",
        "SAT_rass_first_neg_30_last45_nonneg",
    ]

    vent_eligible = (
        cohort_elig[cohort_elig["on_vent_and_sedation"] == 1]
        .sort_values(["hospitalization_id", "recorded_dttm"])
        .reset_index(drop=True)
        .copy()
    )

    for _f in _FLAG_COLS:
        vent_eligible[_f] = np.nan

    # rank_sedation: cumulative zeros per hosp-day
    vent_eligible["rank_sedation"] = np.nan
    for _key, _grp in tqdm(vent_eligible.groupby("hosp_id_day_key"), desc="Rank sedation"):
        _zero_mask = _grp["min_sedation_dose"] == 0
        _ranks = _zero_mask.cumsum() * _zero_mask
        vent_eligible.loc[_grp.index, "rank_sedation"] = _ranks.replace(0, np.nan)

    vent_eligible["rank_sedation_non_ops"] = np.nan
    for _key, _grp in tqdm(vent_eligible.groupby("hosp_id_day_key"), desc="Rank sedation (non-opioid)"):
        _zero_mask = _grp["min_sedation_dose_non_ops"] == 0
        _ranks = _zero_mask.cumsum() * _zero_mask
        vent_eligible.loc[_grp.index, "rank_sedation_non_ops"] = _ranks.replace(0, np.nan)

    vent_eligible["rass"] = vent_eligible["rass"].astype(float)

    # --- Main 6-flag evaluation loop ---
    _delta30 = pd.Timedelta(minutes=30)
    _delta45 = pd.Timedelta(minutes=45)

    for _key, _grp in tqdm(vent_eligible.groupby("hosp_id_day_key"), desc="Evaluating SAT flags (6 methods)"):
        _grp = _grp.sort_values("recorded_dttm")
        _idxs = _grp.index
        _times = _grp["recorded_dttm"].values
        _ranks = _grp["rank_sedation"].values

        for _idx, _cur_time, _rank in zip(_idxs, _times, _ranks):
            if pd.isna(_rank):
                continue

            _fw30 = _grp[(_grp["recorded_dttm"] >= _cur_time) & (_grp["recorded_dttm"] <= _cur_time + _delta30)]
            _fw45 = _grp[(_grp["recorded_dttm"] >= _cur_time) & (_grp["recorded_dttm"] <= _cur_time + _delta45)]
            _pr30 = _grp[(_grp["recorded_dttm"] >= _cur_time - _delta30) & (_grp["recorded_dttm"] < _cur_time)]

            # Check IMV maintained in forward 30 min
            _imv_ok = (_fw30["device_category"].str.lower() == "imv").all() if not _fw30.empty else False
            if not _imv_ok:
                continue

            _flags_set = {}

            # Flag 1: SAT_EHR_delivery — no meds in next 30 min
            _meds_ok = (_fw30[_MED_COLS].isna() | (_fw30[_MED_COLS] == 0)).all().all()
            if _meds_ok:
                _flags_set["SAT_EHR_delivery"] = 1

            # Flag 2: SAT_modified_delivery — no non-opioid meds in next 30 min
            _meds2_ok = (_fw30[_MED_COLS2].isna() | (_fw30[_MED_COLS2] == 0)).all().all()
            if _meds2_ok:
                _flags_set["SAT_modified_delivery"] = 1

            _rass30 = _fw30["rass"].dropna()
            _rass45 = _fw45["rass"].dropna()
            _rass45_ok = not _rass45.empty and _rass45.iloc[-1] >= 0
            _rass30_pre = _pr30["rass"].dropna()

            # Flag 3: SAT_rass_nonneg_30 — all RASS >= 0 in next 30 min
            if not _rass30.empty and (_rass30 >= 0).all():
                _flags_set["SAT_rass_nonneg_30"] = 1

            # Flag 4: SAT_med_halved_rass_pos — all meds <= 50% of prior 30-min max AND last RASS in 45 min >= 0
            if not _pr30.empty and not _fw30.empty:
                _half_max = _pr30[_MED_COLS].max() * 0.5
                _halved_ok = True
                for _med in _MED_COLS:
                    _vals = _fw30[_med].dropna()
                    _vals = _vals[_vals != 0]
                    if not _vals.empty and not (_vals <= _half_max[_med]).all():
                        _halved_ok = False
                        break
                if _halved_ok and _rass45_ok:
                    _flags_set["SAT_med_halved_rass_pos"] = 1

            # Flag 5: SAT_no_meds_rass_pos_45 — no meds in 30 min AND last RASS in 45 min >= 0
            if _meds_ok and _rass45_ok:
                _flags_set["SAT_no_meds_rass_pos_45"] = 1

            # Flag 6: SAT_rass_first_neg_30_last45_nonneg — first prior RASS < 0, last RASS in 45 min >= 0
            if not _rass30_pre.empty and not _rass45.empty:
                if _rass30_pre.iloc[0] < 0 and _rass45.iloc[-1] >= 0:
                    _flags_set["SAT_rass_first_neg_30_last45_nonneg"] = 1

            for _f, _val in _flags_set.items():
                vent_eligible.at[_idx, _f] = _val

    sat_flag_cols = _FLAG_COLS
    return sat_flag_cols, vent_eligible


@app.cell
def _(cohort_elig, np, pl, sat_flag_cols, vent_eligible):
    # Patch: if site has no RASS, zero out RASS-based flags
    if cohort_elig["rass"].nunique() <= 1:
        for _f in ["SAT_rass_nonneg_30", "SAT_med_halved_rass_pos",
                    "SAT_no_meds_rass_pos_45", "SAT_rass_first_neg_30_last45_nonneg"]:
            vent_eligible[_f] = 0
        print("Site has no RASS — RASS-based flags set to 0.")

    cohort_with_delivery = vent_eligible.copy()

    # Day-level flag aggregation
    _max_cols = [
        "sat_screen_pass_fail", "sat_delivery_pass_fail",
        "eligible_event",
    ] + sat_flag_cols
    _agg_pl = (
        pl.from_pandas(cohort_with_delivery[["hosp_id_day_key"] + _max_cols])
        .group_by("hosp_id_day_key")
        .agg([pl.col(c).max().alias(c) for c in _max_cols])
        .sort("hosp_id_day_key")
    )
    df_grouped = _agg_pl.to_pandas()

    df_grouped["sat_flowsheet_delivery_flag"] = np.where(
        ((df_grouped["sat_screen_pass_fail"] == 1) | (df_grouped["sat_delivery_pass_fail"] == 1))
        & (df_grouped["eligible_event"] == 1),
        1, np.nan,
    )

    print(f"Day-level aggregated: {len(df_grouped):,} hosp-days")
    for _col in sat_flag_cols + ["sat_flowsheet_delivery_flag"]:
        _n = (df_grouped[_col] == 1).sum()
        print(f"  {_col}: {_n:,} days flagged")
    return cohort_with_delivery, df_grouped


@app.cell
def _(OUTPUT_SAT, cohort_elig, df_grouped, np, pd, sat_flag_cols):
    from sklearn.metrics import cohen_kappa_score, confusion_matrix
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams["font.family"] = "Arial"
    rcParams["font.size"] = 8

    # Bootstrap kappa CI
    def _kappa_ci_bootstrap(y_true, y_pred, n_boot=2000, ci=0.95, seed=42):
        rng = np.random.default_rng(seed)
        n = len(y_true)
        kappas = []
        for _ in range(n_boot):
            idx = rng.choice(n, size=n, replace=True)
            try:
                kappas.append(cohen_kappa_score(y_true.iloc[idx], y_pred.iloc[idx]))
            except Exception:
                continue
        alpha = (1 - ci) / 2
        return float(np.percentile(kappas, alpha * 100)), float(np.percentile(kappas, (1 - alpha) * 100))

    def _landis_koch(kappa):
        if kappa < 0: return "Poor"
        elif kappa < 0.21: return "Slight"
        elif kappa < 0.41: return "Fair"
        elif kappa < 0.61: return "Moderate"
        elif kappa < 0.81: return "Substantial"
        return "Almost Perfect"

    _con_df = df_grouped.copy()
    _fill = sat_flag_cols + ["sat_flowsheet_delivery_flag"]
    for _c in _fill:
        _con_df[_c] = _con_df[_c].fillna(0)

    _has_rass = cohort_elig["rass"].nunique() > 1
    _metrics_list = []

    for _col in sat_flag_cols:
        if "rass" in _col and not _has_rass:
            continue
        _y_true = _con_df["sat_flowsheet_delivery_flag"]
        _y_pred = _con_df[_col]
        _cm = confusion_matrix(_y_true, _y_pred)
        _tn, _fp, _fn, _tp = _cm.ravel()
        _total = _cm.sum()
        _accuracy = (_tp + _tn) / _total
        _precision = _tp / (_tp + _fp) if _tp + _fp else 0
        _recall = _tp / (_tp + _fn) if _tp + _fn else 0
        _f1 = 2 * _precision * _recall / (_precision + _recall) if _precision + _recall else 0
        _specificity = _tn / (_tn + _fp) if _tn + _fp else 0
        _kappa = cohen_kappa_score(_y_true, _y_pred)
        _kappa_lo, _kappa_hi = _kappa_ci_bootstrap(_y_true, _y_pred)

        # JAMA-style confusion matrix plot
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
        _fig.savefig(OUTPUT_SAT / f"confusion_matrix_{_col}.png", bbox_inches="tight", dpi=300)
        plt.close(_fig)

        _metrics_list.append({
            "Column": _col, "TP": _tp, "FP": _fp, "FN": _fn, "TN": _tn,
            "Accuracy": _accuracy, "Precision": _precision, "Recall": _recall,
            "F1": _f1, "Specificity": _specificity,
            "Cohen_Kappa": _kappa, "Kappa_CI_lower": _kappa_lo, "Kappa_CI_upper": _kappa_hi,
            "Kappa_Interpretation": _landis_koch(_kappa),
        })

    concordance_df = pd.DataFrame(_metrics_list)
    concordance_df.to_csv(OUTPUT_SAT / "delivery_concordance_summary.csv", index=False)
    print(concordance_df.to_string(index=False))
    return


@app.cell
def _(
    OUTPUT_SAT,
    SITE_NAME,
    cohort,
    cohort_with_delivery,
    df_grouped,
    pd,
    sat_flag_cols,
):
    # Get hospital_id from cohort
    _hosp_ids = cohort.select("hospitalization_id", "hospital_id").unique().to_pandas() if "hospital_id" in cohort.columns else None

    if _hosp_ids is not None:
        _final = df_grouped.copy()
        _hid = cohort_with_delivery[["hospitalization_id", "hosp_id_day_key"]].drop_duplicates()
        _final = _final.merge(_hid, on="hosp_id_day_key", how="left")
        _final = _final.merge(_hosp_ids, on="hospitalization_id", how="left")

        _rows = []
        for _hosp in _final["hospital_id"].dropna().unique():
            _sub = _final[_final["hospital_id"] == _hosp]
            _eligible_n = int((_sub["eligible_event"] == 1).sum())
            _row = {"Site_Hospital": f"{SITE_NAME}_{_hosp}", "eligible_event_count": _eligible_n}
            for _flag in sat_flag_cols + ["sat_flowsheet_delivery_flag"]:
                _n = int((_sub[_flag] == 1).sum())
                _row[f"{_flag}_count"] = _n
                _row[f"pct_{_flag}"] = round(_n / _eligible_n * 100, 2) if _eligible_n > 0 else 0.0
            _rows.append(_row)
        hospital_summary = pd.DataFrame(_rows)
    else:
        # Single-site: summarize overall
        _eligible_n = int((df_grouped["eligible_event"] == 1).sum())
        _row = {"Site_Hospital": SITE_NAME, "eligible_event_count": _eligible_n}
        for _flag in sat_flag_cols + ["sat_flowsheet_delivery_flag"]:
            _n = int((df_grouped[_flag] == 1).sum())
            _row[f"{_flag}_count"] = _n
            _row[f"pct_{_flag}"] = round(_n / _eligible_n * 100, 2) if _eligible_n > 0 else 0.0
        hospital_summary = pd.DataFrame([_row])

    hospital_summary.to_csv(OUTPUT_SAT / f"sat_stats_{SITE_NAME}.csv", index=False)
    print(hospital_summary.T)
    return


@app.cell
def _(OUTPUT_SAT, df_grouped, pl):
    _out = pl.from_pandas(df_grouped)
    _path = OUTPUT_SAT / "sat_day_level.parquet"
    _out.write_parquet(str(_path))
    print(f"Saved: {_path} ({_out.height:,} rows x {_out.width} cols)")
    return


@app.cell
def _(df_grouped):
    df_grouped
    return


if __name__ == "__main__":
    app.run()
