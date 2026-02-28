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
    OUTPUT_SBT = OUTPUT_PHI / "sbt_standard"
    OUTPUT_SBT.mkdir(parents=True, exist_ok=True)

    # SBT constants (from definitions_source_of_truth.py)
    VENT_DAY_ANCHOR_HOUR = 6
    SBT_CONTROLLED_MODES = [
        "assist control-volume control",
        "pressure control",
        "pressure-regulated volume control",
        "simv",
    ]
    SBT_MIN_CONTROLLED_MODE_HOURS = 12
    SBT_PS_MAX = 8
    SBT_PEEP_MAX = 8
    SBT_PRIMARY_DURATION_MIN = 2
    SBT_SECONDARY_DURATIONS = [5, 30]
    SBT_FIO2_MAX = 0.50
    SBT_SPO2_MIN = 88
    SBT_NEE_MAX = 0.2

    print(f"Site: {SITE_NAME}")
    print(f"Output: {OUTPUT_SBT}")
    return (
        OUTPUT_PHI, OUTPUT_SBT, SITE_NAME,
        SBT_CONTROLLED_MODES, SBT_FIO2_MAX, SBT_MIN_CONTROLLED_MODE_HOURS,
        SBT_NEE_MAX, SBT_PEEP_MAX, SBT_PRIMARY_DURATION_MIN, SBT_PS_MAX,
        SBT_SECONDARY_DURATIONS, SBT_SPO2_MIN, VENT_DAY_ANCHOR_HOUR,
        np, pd, pl, tqdm,
    )


@app.cell
def _(OUTPUT_PHI, VENT_DAY_ANCHOR_HOUR, pd, pl):
    # Load wide dataset and cohort
    wide = pl.read_parquet(OUTPUT_PHI / "wide_dataset.parquet")
    print(f"Wide dataset: {wide.height:,} rows x {wide.width} cols")

    cohort = pl.read_parquet(OUTPUT_PHI / "cohort.parquet")
    print(f"Cohort: {cohort.height:,} hospitalizations")

    # Convert to pandas for group-based sequential processing
    df_raw = wide.to_pandas()
    df_raw["recorded_dttm"] = pd.to_datetime(df_raw["recorded_dttm"])
    df_raw["first_icu_start"] = pd.to_datetime(df_raw["first_icu_start"])
    df_raw["first_icu_end"] = pd.to_datetime(df_raw["first_icu_end"])

    # 06:00-anchored vent-day date
    df_raw["vent_day_date"] = (
        df_raw["recorded_dttm"] - pd.Timedelta(hours=VENT_DAY_ANCHOR_HOUR)
    ).dt.normalize()
    df_raw["hosp_id_day_key"] = (
        df_raw["hospitalization_id"].astype(str) + "_" + df_raw["vent_day_date"].dt.strftime("%Y-%m-%d")
    )

    # Join extubation_time from cohort
    _cohort_pd = cohort.to_pandas()
    _cohort_pd["extubation_time"] = pd.to_datetime(_cohort_pd["extubation_time"])
    _ext = _cohort_pd[["hospitalization_id", "extubation_time"]].drop_duplicates()
    df_raw = df_raw.merge(_ext, on="hospitalization_id", how="left")

    # extubated flag: 1 if recorded_dttm >= extubation_time
    df_raw["extubated"] = (
        (df_raw["extubation_time"].notna()) & (df_raw["recorded_dttm"] >= df_raw["extubation_time"])
    ).astype(int)

    df_raw = df_raw.sort_values(["hospitalization_id", "recorded_dttm"]).reset_index(drop=True)

    print(f"Preprocessed: {len(df_raw):,} rows, {df_raw['hosp_id_day_key'].nunique():,} hosp-day keys")
    print(f"Extubation times joined: {df_raw['extubation_time'].notna().sum():,} rows have extubation_time")
    return (df_raw,)


@app.cell
def _(
    SBT_CONTROLLED_MODES, SBT_FIO2_MAX, SBT_NEE_MAX, SBT_PEEP_MAX,
    SBT_SPO2_MIN, df_raw,
):
    # Stability flags + vent day classification (vectorized)
    df_classified = df_raw.copy()
    df_classified["device_category"] = df_classified["device_category"].astype(str).str.lower()
    df_classified["mode_category"] = df_classified["mode_category"].astype(str).str.lower()

    # IMV flag: device_category == 'imv' AND mode_category in controlled modes
    df_classified["imv_flag"] = (
        (df_classified["device_category"] == "imv")
        & (df_classified["mode_category"].isin(SBT_CONTROLLED_MODES))
    )

    # Hemodynamic stability: NEE <= 0.2
    df_classified["hemodynamic_stability"] = (df_classified["nee"].fillna(0) <= SBT_NEE_MAX).astype(int)

    # Respiratory stability: FiO2 <= 0.5 AND PEEP <= 8 AND SpO2 >= 88
    df_classified["respiratory_stability"] = (
        (df_classified["fio2_set"].fillna(1.0) <= SBT_FIO2_MAX)
        & (df_classified["peep_set"].fillna(99) <= SBT_PEEP_MAX)
        & (df_classified["spo2"].fillna(0) >= SBT_SPO2_MIN)
    ).astype(int)

    # Vent day: any row in that hosp_id_day_key has device_category == 'imv'
    _imv_days = df_classified[df_classified["device_category"] == "imv"]["hosp_id_day_key"].unique()
    df_classified["vent_day"] = df_classified["hosp_id_day_key"].isin(_imv_days).astype(int)

    # Vent day without paralytics
    _para_days = (
        df_classified[df_classified["vent_day"] == 1]
        .groupby("hosp_id_day_key")["max_paralytics"]
        .max()
    )
    _no_para_days = set(_para_days[_para_days == 0].index)
    df_classified["vent_day_without_paralytics"] = (
        (df_classified["vent_day"] == 1) & df_classified["hosp_id_day_key"].isin(_no_para_days)
    ).astype(int)

    print("Stability flags computed:")
    print(f"  imv_flag=True: {df_classified['imv_flag'].sum():,}")
    print(f"  hemodynamic_stability=1: {df_classified['hemodynamic_stability'].sum():,}")
    print(f"  respiratory_stability=1: {df_classified['respiratory_stability'].sum():,}")
    _vent_days = df_classified[df_classified["vent_day"] == 1]["hosp_id_day_key"].nunique()
    _vent_no_para = df_classified[df_classified["vent_day_without_paralytics"] == 1]["hosp_id_day_key"].nunique()
    print(f"Vent days (at least 1 IMV): {_vent_days:,}")
    print(f"Vent days w/o paralytics: {_vent_no_para:,}")
    return (df_classified,)


@app.cell
def _(SBT_MIN_CONTROLLED_MODE_HOURS, VENT_DAY_ANCHOR_HOUR, df_classified, pd, tqdm):
    # Condition 1: 12h contiguous controlled-mode IMV check
    df_eligible = df_classified.copy()

    _cond1_threshold = pd.Timedelta(hours=SBT_MIN_CONTROLLED_MODE_HOURS)
    _cond1_window_start_offset = pd.Timedelta(hours=VENT_DAY_ANCHOR_HOUR) - pd.Timedelta(days=1)
    _cond1_window_end_offset = pd.Timedelta(hours=VENT_DAY_ANCHOR_HOUR)

    df_eligible["IMV_Controlled_met_time"] = pd.NaT
    df_eligible["eligible_day"] = 0

    # Build hospitalization lookup for 24h lookback
    _hosp_groups = {
        hosp_id: grp.copy().sort_values("recorded_dttm")
        for hosp_id, grp in df_eligible.groupby("hospitalization_id")
    }

    # Only process vent days without paralytics
    _candidates = df_eligible[df_eligible["vent_day_without_paralytics"] == 1]
    _groups = _candidates.groupby(["hospitalization_id", "vent_day_date"])

    for (hosp_id, curr_day), day_group in tqdm(_groups, desc="SBT eligibility (12h controlled mode)"):
        cond1_start = curr_day + _cond1_window_start_offset
        cond1_end = curr_day + _cond1_window_end_offset

        hosp_df = _hosp_groups[hosp_id]
        cond1_df = hosp_df[
            (hosp_df["recorded_dttm"] >= cond1_start)
            & (hosp_df["recorded_dttm"] <= cond1_end)
        ].copy()

        if cond1_df.empty or not cond1_df["imv_flag"].any():
            continue

        # Check paralytics in the 24h lookback window
        if cond1_df["max_paralytics"].max() > 0:
            continue

        # Find contiguous segments where imv_flag is True
        cond1_df["seg"] = (cond1_df["imv_flag"] != cond1_df["imv_flag"].shift()).cumsum()
        valid_segs = cond1_df[cond1_df["imv_flag"]].groupby("seg")

        for seg_id, seg_df in valid_segs:
            seg_df = seg_df.sort_values("recorded_dttm").copy()
            seg_df["duration"] = seg_df["recorded_dttm"].diff().fillna(pd.Timedelta(seconds=0))
            seg_df["cum_duration"] = seg_df["duration"].cumsum()
            if seg_df["cum_duration"].iloc[-1] >= _cond1_threshold:
                flag_row = seg_df[seg_df["cum_duration"] >= _cond1_threshold].iloc[0]
                flag_idx = flag_row.name
                flag_time = flag_row["recorded_dttm"]
                df_eligible.loc[flag_idx, "IMV_Controlled_met_time"] = flag_time
                df_eligible.loc[day_group.index, "eligible_day"] = 1
                break

    _eligible_days = df_eligible[df_eligible["eligible_day"] == 1]["hosp_id_day_key"].nunique()
    _vent_no_para = df_eligible[df_eligible["vent_day_without_paralytics"] == 1]["hosp_id_day_key"].nunique()
    _pct = (_eligible_days / _vent_no_para * 100) if _vent_no_para > 0 else 0
    print(f"Eligible days: {_eligible_days:,} / {_vent_no_para:,} ({_pct:.1f}%)")
    return (df_eligible,)


@app.cell
def _(SBT_PEEP_MAX, SBT_PRIMARY_DURATION_MIN, SBT_PS_MAX, SBT_SECONDARY_DURATIONS, df_eligible, np, pd):
    # Flip detection (SBT delivery)
    df_delivery = df_eligible.copy()
    _durations_min = [SBT_PRIMARY_DURATION_MIN] + SBT_SECONDARY_DURATIONS  # [2, 5, 30]

    # Compute flip_check_flag vectorized
    _mode_cat = df_delivery["mode_category"].fillna("")
    _mode_name = df_delivery["mode_name"].astype(str).str.lower() if "mode_name" in df_delivery.columns else pd.Series("", index=df_delivery.index)
    _cond_imv = df_delivery["device_category"] == "imv"
    _cond_mode_ps = _mode_cat.str.contains("pressure support|cpap", regex=True, na=False)
    _cond_ps_le = df_delivery["pressure_support_set"].fillna(99) <= SBT_PS_MAX
    _cond_peep_le = df_delivery["peep_set"].fillna(99) <= SBT_PEEP_MAX
    _conditionA = _cond_mode_ps & _cond_ps_le & _cond_peep_le
    _cond_tpiece = _mode_name.str.match(r"^t[-\s]?piece$", na=False)
    _composite = _conditionA | _cond_tpiece
    _passed = _cond_imv & _composite

    df_delivery["flip_check_flag"] = False
    _mask_eligible = df_delivery["eligible_day"] == 1
    df_delivery.loc[_mask_eligible, "flip_check_flag"] = _passed[_mask_eligible]

    # Initialize delivery and diagnostic columns
    for _d in _durations_min:
        df_delivery[f"EHR_Delivery_{_d}mins"] = np.nan
    df_delivery["first_flip_time"] = pd.NaT
    df_delivery["flip_skip_reason"] = None

    # Compute min IMV_Controlled_met_time per eligible group
    df_delivery.loc[_mask_eligible, "min_met_time"] = (
        df_delivery.loc[_mask_eligible]
        .groupby(["hospitalization_id", "vent_day_date"])["IMV_Controlled_met_time"]
        .transform("min")
    )

    # Per-group flip processing (ported from pySBT.process_diagnostic_flip_sbt_optimized_v2)
    def _process_flip_group(group):
        group = group.sort_values("recorded_dttm").copy()
        n = len(group)
        if n == 0:
            return group

        times = group["recorded_dttm"].values.astype("datetime64[ns]")
        flip_int = group["flip_check_flag"].astype(int).values

        def compute_sustained(delta_minutes):
            delta = np.timedelta64(delta_minutes, "m")
            boundaries = np.searchsorted(times, times + delta, side="right")
            cnt_total = boundaries - np.arange(n)
            cumsum = np.cumsum(flip_int)
            cnt_pass = np.empty(n, dtype=int)
            for i in range(n):
                end = boundaries[i] - 1
                if end < i:
                    cnt_pass[i] = 0
                else:
                    cnt_pass[i] = cumsum[end] - (cumsum[i - 1] if i > 0 else 0)
            return (cnt_total == cnt_pass) & group["flip_check_flag"].values

        sustained = {}
        for d in _durations_min:
            sustained[d] = compute_sustained(d)

        # Primary duration logic
        primary_duration = _durations_min[0]
        candidate_indices = group.index[group["flip_check_flag"]].tolist()
        for idx in candidate_indices:
            row_pos = group.index.get_loc(idx)
            group.at[idx, "first_flip_time"] = group.at[idx, "recorded_dttm"]
            if group.at[idx, "recorded_dttm"] < group.at[idx, "min_met_time"]:
                group.at[idx, "flip_skip_reason"] = "Flip before IMV_Controlled_met_time"
                continue
            if sustained[primary_duration][row_pos]:
                group.at[idx, f"EHR_Delivery_{primary_duration}mins"] = 1
                group.at[idx, "flip_skip_reason"] = None
                break
            else:
                group.at[idx, "flip_skip_reason"] = f"ehr_delivery_{primary_duration}min not possible"

        # Secondary durations (independently)
        for d in _durations_min[1:]:
            for idx in candidate_indices:
                row_pos = group.index.get_loc(idx)
                if group.at[idx, "recorded_dttm"] < group.at[idx, "min_met_time"]:
                    continue
                if sustained[d][row_pos]:
                    group.at[idx, f"EHR_Delivery_{d}mins"] = 1
                    break

        return group

    _eligible_df = df_delivery[_mask_eligible].copy()
    _processed = (
        _eligible_df
        .groupby(["hospitalization_id", "vent_day_date"], group_keys=False)
        .apply(_process_flip_group)
    )
    df_delivery.update(_processed)

    # Clean up helper columns
    df_delivery.drop(columns=["min_met_time"], inplace=True, errors="ignore")

    _n_2min = (df_delivery["EHR_Delivery_2mins"] == 1).sum()
    _n_5min = (df_delivery["EHR_Delivery_5mins"] == 1).sum()
    _n_30min = (df_delivery["EHR_Delivery_30mins"] == 1).sum()
    print(f"EHR Delivery 2min: {_n_2min:,}")
    print(f"EHR Delivery 5min: {_n_5min:,}")
    print(f"EHR Delivery 30min: {_n_30min:,}")
    return (df_delivery,)


@app.cell
def _(df_delivery, np, pd):
    # Post-delivery: extubation tracking
    df_annotated = df_delivery.copy()
    df_annotated["flag_2_45_extubated"] = np.nan
    df_annotated["delta_to_extubation_mins"] = np.nan

    _group_cols = ["hospitalization_id", "vent_day_date"]
    for (hosp_id, day), group in df_annotated.groupby(_group_cols):
        group = group.sort_values("recorded_dttm")

        # flag_2_45_extubated: extubation within 45 min of 2-min flip
        flip_row = group[
            (group["EHR_Delivery_2mins"] == 1) & (group["first_flip_time"].notna())
        ]
        if not flip_row.empty:
            flip_time = flip_row.iloc[0]["first_flip_time"]
            window_end = flip_time + pd.Timedelta(minutes=45)
            ext_mask = (
                (group["recorded_dttm"] > flip_time)
                & (group["recorded_dttm"] <= window_end)
                & (group["extubated"] == 1)
            )
            if ext_mask.any():
                df_annotated.loc[flip_row.index[0], "flag_2_45_extubated"] = 1

        # delta_to_extubation_mins: time from 30-min flip to first extubation
        flip_30 = group[
            (group["EHR_Delivery_30mins"] == 1) & (group["first_flip_time"].notna())
        ]
        if not flip_30.empty:
            flip_time_30 = flip_30.iloc[0]["first_flip_time"]
            post_ext = group[
                (group["recorded_dttm"] > flip_time_30) & (group["extubated"] == 1)
            ]
            if not post_ext.empty:
                ext_time = post_ext.iloc[0]["recorded_dttm"]
                delta = (ext_time - flip_time_30).total_seconds() / 60.0
                df_annotated.loc[flip_30.index[0], "delta_to_extubation_mins"] = delta

    _n_ext45 = (df_annotated["flag_2_45_extubated"] == 1).sum()
    _n_delta = df_annotated["delta_to_extubation_mins"].notna().sum()
    print(f"Extubation within 45min of 2-min flip: {_n_ext45:,}")
    print(f"Delta to extubation computed: {_n_delta:,}")
    return (df_annotated,)


@app.cell
def _(df_annotated, np):
    # Day-level aggregation & summary stats
    # Convert EHR delivery datetime columns to binary int for aggregation
    _df = df_annotated.copy()
    for _col in ["EHR_Delivery_2mins", "EHR_Delivery_5mins", "EHR_Delivery_30mins"]:
        _df[_col] = _df[_col].notna().astype(int)
        _df[_col] = _df[_col].where(_df[_col] == 1, np.nan)

    # Fill forward flip_skip_reason within each hosp-day
    _df["flip_skip_reason"] = _df.groupby("hosp_id_day_key")["flip_skip_reason"].transform(
        lambda x: x.ffill().bfill()
    )

    grouped_df = (
        _df.groupby("hosp_id_day_key")
        .agg(
            hospitalization_id=("hospitalization_id", "first"),
            eligible_day=("eligible_day", "max"),
            vent_day=("vent_day", "max"),
            vent_day_without_paralytics=("vent_day_without_paralytics", "max"),
            EHR_Delivery_2mins=("EHR_Delivery_2mins", "max"),
            EHR_Delivery_5mins=("EHR_Delivery_5mins", "max"),
            EHR_Delivery_30mins=("EHR_Delivery_30mins", "max"),
            sbt_screen_pass_fail=("sbt_screen_pass_fail", "max"),
            sbt_delivery_pass_fail=("sbt_delivery_pass_fail", "max"),
            extubated=("extubated", "max"),
            flag_2_45_extubated=("flag_2_45_extubated", "max"),
            flip_skip_reason=("flip_skip_reason", lambda x: x.dropna().iloc[-1] if x.dropna().size > 0 else np.nan),
        )
        .reset_index()
    )

    _total = len(grouped_df)
    _vent = (grouped_df["vent_day"] == 1).sum()
    _vent_no_para = (grouped_df["vent_day_without_paralytics"] == 1).sum()
    _eligible = (grouped_df["eligible_day"] == 1).sum()
    _ehr_2 = (grouped_df["EHR_Delivery_2mins"] == 1).sum()
    _ehr_5 = (grouped_df["EHR_Delivery_5mins"] == 1).sum()
    _ehr_30 = (grouped_df["EHR_Delivery_30mins"] == 1).sum()
    _sbt_s = (grouped_df["sbt_screen_pass_fail"] == 1).sum()
    _sbt_d = (grouped_df["sbt_delivery_pass_fail"] == 1).sum()
    _ext = (grouped_df["extubated"] == 1).sum()
    _ext45 = (grouped_df["flag_2_45_extubated"] == 1).sum()

    print("=== Day-Level Summary ===")
    print(f"Total hosp-days: {_total:,}")
    print(f"Vent days: {_vent:,}")
    print(f"Vent days w/o paralytics: {_vent_no_para:,}")
    print(f"Eligible days: {_eligible:,}")
    print(f"EHR Delivery 2min: {_ehr_2:,}")
    print(f"EHR Delivery 5min: {_ehr_5:,}")
    print(f"EHR Delivery 30min: {_ehr_30:,}")
    print(f"SBT screen pass (flowsheet): {_sbt_s:,}")
    print(f"SBT delivery pass (flowsheet): {_sbt_d:,}")
    print(f"Extubated: {_ext:,}")
    print(f"Extubated within 45min of flip: {_ext45:,}")
    return (grouped_df,)


@app.cell
def _(OUTPUT_SBT, grouped_df, np, pd):
    # Concordance: EHR delivery columns vs sbt_delivery_pass_fail
    from sklearn.metrics import cohen_kappa_score, confusion_matrix
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams["font.family"] = "Arial"
    rcParams["font.size"] = 8

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

    # Prepare concordance dataframe (eligible days only)
    _con_df = grouped_df[grouped_df["eligible_day"] == 1].copy()

    _ehr_cols = ["EHR_Delivery_2mins", "EHR_Delivery_5mins", "EHR_Delivery_30mins"]
    _ref_col = "sbt_delivery_pass_fail"

    for _c in _ehr_cols + [_ref_col]:
        _con_df[_c] = _con_df[_c].fillna(0).astype(int)

    _metrics_list = []
    for _col in _ehr_cols:
        _y_true = _con_df[_ref_col]
        _y_pred = _con_df[_col]
        _cm = confusion_matrix(_y_true, _y_pred, labels=[0, 1])
        _tn, _fp, _fn, _tp = _cm.ravel()
        _total = _cm.sum()
        _accuracy = (_tp + _tn) / _total if _total else 0
        _precision = _tp / (_tp + _fp) if (_tp + _fp) else 0
        _recall = _tp / (_tp + _fn) if (_tp + _fn) else 0
        _f1 = 2 * _precision * _recall / (_precision + _recall) if (_precision + _recall) else 0
        _specificity = _tn / (_tn + _fp) if (_tn + _fp) else 0
        _kappa = cohen_kappa_score(_y_true, _y_pred)
        _kappa_lo, _kappa_hi = _kappa_ci_bootstrap(_y_true, _y_pred)

        # Confusion matrix plot
        _cm_pct = _cm / _total * 100 if _total else _cm * 0
        _fig, _ax = plt.subplots(figsize=(3.5, 3.0))
        _ax.imshow(_cm, cmap="cividis")
        _ax.set_xticks([0, 1]); _ax.set_yticks([0, 1])
        _ax.set_xticklabels(["No Delivery", "Delivery"])
        _ax.set_yticklabels(["No Delivery", "Delivery"])
        _ax.set_xlabel(f"{_col}"); _ax.set_ylabel("Flowsheet SBT delivery")
        _ax.set_title(f"Concordance: flowsheet vs {_col}", fontweight="bold")
        _ax.spines["top"].set_visible(False); _ax.spines["right"].set_visible(False)
        for _i in range(2):
            for _j in range(2):
                _ax.text(_j, _i, f"{_cm[_i, _j]}\n({_cm_pct[_i, _j]:.1f}%)",
                         ha="center", va="center",
                         color="white" if _cm[_i, _j] > _cm.max() / 2 else "black")
        _fig.tight_layout()
        _fig.savefig(OUTPUT_SBT / f"confusion_matrix_{_col}.png", bbox_inches="tight", dpi=300)
        plt.close(_fig)

        _metrics_list.append({
            "Column": _col, "TP": _tp, "FP": _fp, "FN": _fn, "TN": _tn,
            "Accuracy": _accuracy, "Precision": _precision, "Recall": _recall,
            "F1": _f1, "Specificity": _specificity,
            "Cohen_Kappa": _kappa, "Kappa_CI_lower": _kappa_lo, "Kappa_CI_upper": _kappa_hi,
            "Kappa_Interpretation": _landis_koch(_kappa),
        })

    concordance_df = pd.DataFrame(_metrics_list)
    concordance_df.to_csv(OUTPUT_SBT / "delivery_concordance_summary.csv", index=False)
    print(concordance_df.to_string(index=False))
    return


@app.cell
def _(OUTPUT_SBT, df_annotated, grouped_df, pl):
    # Save results
    _day_level = pl.from_pandas(grouped_df)
    _day_path = OUTPUT_SBT / "sbt_results.parquet"
    _day_level.write_parquet(str(_day_path))
    _day_size = _day_path.stat().st_size / 1024 / 1024
    print(f"Saved day-level: {_day_path} ({_day_level.height:,} rows, {_day_size:.1f} MB)")

    # Also save to output_phi root for downstream use
    _root_path = OUTPUT_SBT.parent / "sbt_results.parquet"
    _day_level.write_parquet(str(_root_path))
    print(f"Saved copy: {_root_path}")

    # Save row-level annotated dataset
    _row_path = OUTPUT_SBT / "sbt_annotated_wide.parquet"
    _row_level = pl.from_pandas(df_annotated)
    _row_level.write_parquet(str(_row_path))
    _row_size = _row_path.stat().st_size / 1024 / 1024
    print(f"Saved row-level: {_row_path} ({_row_level.height:,} rows, {_row_size:.1f} MB)")
    return


@app.cell
def _(grouped_df):
    grouped_df
    return


if __name__ == "__main__":
    app.run()
