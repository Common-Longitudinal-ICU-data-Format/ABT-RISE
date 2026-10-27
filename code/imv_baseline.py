"""Replay the current event cell and run a byte-identical, privately staged notebook.

No detection expressions in code/01_cohort.py are rewritten. The staging boundary
changes file locations, not clinical criteria. Patient-level outputs stay local.
"""
from __future__ import annotations

import argparse
import ast
from contextlib import redirect_stderr, redirect_stdout
import gc
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any

import polars as pl

EVENT_FIELDS = [
    "hospitalization_id", "intubation_time", "extubation_time",
    "extubation_method", "no_extubation_reason", "imv_duration_hours",
]
RESP_FIELDS = [
    "hospitalization_id", "recorded_dttm", "device_category", "mode_category",
    "tracheostomy", "lpm_set", "peep_set", "tidal_volume_set", "resp_rate_set",
    "pressure_control_set", "pressure_support_set",
]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def private_dir(path: Path) -> Path:
    path = Path(path)
    old_umask = os.umask(0o077)
    try:
        path.mkdir(parents=True, exist_ok=True)
    finally:
        os.umask(old_umask)
    if os.name == "posix":
        path.chmod(0o700)
    return path


def private_json(path: Path, value: Any) -> None:
    private_dir(path.parent)
    path.write_text(json.dumps(value, indent=2, default=str) + "\n")
    if os.name == "posix":
        path.chmod(0o600)


def private_parquet(path: Path, frame: pl.DataFrame) -> None:
    private_dir(path.parent)
    frame.write_parquet(path)
    if os.name == "posix":
        path.chmod(0o600)


def replay_baseline(
    respiratory: pl.DataFrame,
    hospitals: pl.DataFrame,
    source: Path | None = None,
    *,
    expected_sha256: str | None = None,
    return_branch_candidates: bool = False,
) -> pl.DataFrame | tuple[pl.DataFrame, pl.DataFrame]:
    """Execute unchanged source; optionally capture proven same-time branch ties."""
    source = Path(source or Path(__file__).with_name("01_cohort.py"))
    if expected_sha256 is not None and file_sha256(source) != expected_sha256:
        raise ValueError("Baseline source hash differs from the requested snapshot")
    tree = ast.parse(source.read_text(), filename=str(source))
    matches = [
        node for node in tree.body if isinstance(node, ast.FunctionDef)
        and [arg.arg for arg in node.args.args] == ["hosp_final", "pl", "resp_wf"]
    ]
    if len(matches) != 1:
        raise ValueError("Expected exactly one original baseline event cell")
    # Removing its registration decorator does not change its statements/expressions.
    cell = matches[0]
    cell.decorator_list = []
    module = ast.Module(body=[cell], type_ignores=[])
    namespace: dict[str, Any] = {}
    exec(compile(ast.fix_missing_locations(module), str(source), "exec"), namespace)
    function = namespace[cell.name]
    captured: dict[str, pl.DataFrame] = {}

    def capture_return(frame: Any, event: str, arg: Any) -> None:
        if event == "return" and frame.f_code is function.__code__:
            for name in ("_first_extubation", "_extub_lpm"):
                if name in frame.f_locals:
                    captured[name] = frame.f_locals[name]

    previous_profile = sys.getprofile()
    try:
        if return_branch_candidates:
            sys.setprofile(capture_return)
        with redirect_stdout(io.StringIO()):
            result = function(hospitals, pl, respiratory)
    finally:
        if return_branch_candidates:
            sys.setprofile(previous_profile)
    if not return_branch_candidates:
        return result[0]
    if set(captured) != {"_first_extubation", "_extub_lpm"}:
        raise ValueError("Cannot capture the original baseline branch candidates")
    columns = ["hospitalization_id", "extubation_time"]
    ties = (
        captured["_first_extubation"].select(columns)
        .join(captured["_extub_lpm"].select(columns), on="hospitalization_id", how="inner", suffix="_flow")
        .filter(pl.col("extubation_time").is_not_null() & (pl.col("extubation_time") == pl.col("extubation_time_flow")))
        .select(columns).unique()
    )
    return result[0], ties


def assert_baseline_parity(
    reference: pl.DataFrame, replay: pl.DataFrame, *, method_ties: pl.DataFrame | None = None,
) -> dict[str, int]:
    """Strict event-time parity; only independently captured same-time labels may vary."""
    missing = set(EVENT_FIELDS) - set(reference.columns) | set(EVENT_FIELDS) - set(replay.columns)
    if missing:
        raise ValueError("Baseline parity input is missing required event columns")
    keys = reference.select("hospitalization_id")
    if keys.n_unique() != keys.height or replay["hospitalization_id"].n_unique() != replay.height:
        raise ValueError("Baseline parity requires unique hospitalization keys")
    if replay.join(keys, on="hospitalization_id", how="anti").height:
        raise RuntimeError("Baseline replay contains unexpected hospitalization keys")
    actual = keys.join(replay.select(EVENT_FIELDS), on="hospitalization_id", how="left")
    expected = reference.select(EVENT_FIELDS)
    expected = expected.sort("hospitalization_id")
    actual = actual.sort("hospitalization_id")
    core_fields = [name for name in EVENT_FIELDS if name != "extubation_method"]
    if not expected.select(core_fields).equals(actual.select(core_fields)):
        raise RuntimeError("Baseline event replay differs from the original notebook output")
    paired = expected.join(actual.select("hospitalization_id", "extubation_method"), on="hospitalization_id", suffix="_replay")
    differing_labels = paired.filter(~pl.col("extubation_method").eq_missing(pl.col("extubation_method_replay")))
    if differing_labels.height:
        allowed = ["2-row IMV", "2-row IMV+LPM"]
        if method_ties is None or differing_labels.filter(
            ~pl.col("extubation_method").is_in(allowed) | ~pl.col("extubation_method_replay").is_in(allowed)
        ).height or differing_labels.join(method_ties, on=["hospitalization_id", "extubation_time"], how="anti").height:
            raise RuntimeError("Baseline method labels differ without proven equal branch timestamps")
    return {"rows_compared": expected.height, "mismatches": 0, "method_tie_disagreements": differing_labels.height}


def _snapshot_originals(project: Path, excluded: Path) -> dict[str, str]:
    paths = [project / "code" / "01_cohort.py", project / "code" / "helper.py", project / "clif_config.json"]
    for name in ("output", "output_phi", "output_to_share"):
        directory = project / name
        if directory.exists():
            paths.extend(p for p in directory.rglob("*") if p.is_file() and not p.resolve().is_relative_to(excluded.resolve()))
    return {str(p.relative_to(project)): file_sha256(p) for p in paths if p.is_file()}


def stage_baseline(project: Path, output_root: Path, *, waterfall_cache: Path | None = None) -> Path:
    """Create an isolated project using identical source and an independent cache copy."""
    project = Path(project).resolve()
    output_root = private_dir(Path(output_root).resolve())
    staged = Path(tempfile.mkdtemp(prefix="baseline-", dir=output_root))
    private_dir(staged / "code")
    for filename in ("01_cohort.py", "helper.py"):
        source = project / "code" / filename
        destination = staged / "code" / filename
        shutil.copyfile(source, destination)
        destination.chmod(0o600)
    config = json.loads((project / "clif_config.json").read_text())
    data_directory = Path(config["data_directory"]).expanduser()
    if not data_directory.is_absolute():
        data_directory = project / data_directory
    config["data_directory"] = str(data_directory.resolve())
    private_json(staged / "clif_config.json", config)
    cache = Path(waterfall_cache) if waterfall_cache is not None else project / "output_phi" / "resp_waterfall_cohort.parquet"
    if waterfall_cache is not None and not cache.is_file():
        raise ValueError("Requested baseline waterfall cache does not exist")
    if cache.is_file():
        private_dir(staged / "output_phi")
        destination = staged / "output_phi" / "resp_waterfall_cohort.parquet"
        shutil.copyfile(cache, destination)
        destination.chmod(0o600)
    return staged


def _eligible_panel(definitions: dict[str, Any]) -> pl.DataFrame:
    """Retain the incumbent's pre-IMV eligibility, including its early airway exclusions."""
    hospitals = definitions["hosp_24h"]
    respiratory = definitions["rs_all_df"]
    early = (
        respiratory.join(hospitals.select("hospitalization_id", "admission_dttm"), on="hospitalization_id", how="inner")
        .filter(pl.col("recorded_dttm").is_between(pl.col("admission_dttm"), pl.col("admission_dttm") + pl.duration(hours=24)))
        .filter((pl.col("tracheostomy") == 1) | (pl.col("device_category").str.to_lowercase() == "trach collar"))
        .select("hospitalization_id").unique()
    )
    return (
        hospitals.join(early, on="hospitalization_id", how="anti")
        .join(definitions["icu_24h"], on="hospitalization_id", how="inner", suffix="_icu")
        .select("hospitalization_id", "patient_id", "admission_dttm", "discharge_dttm", "age_at_admission", "discharge_category", "icu_in_dttm", "icu_out_dttm")
        .unique(subset=["hospitalization_id", "icu_in_dttm", "icu_out_dttm"])
        .sort("hospitalization_id", "icu_in_dttm")
    )


def run_baseline(project: Path, output_root: Path, *, waterfall_cache: Path | None = None) -> dict[str, Any]:
    """Run the original app privately, capture original observations, and prove parity."""
    project = Path(project).resolve()
    output_root = private_dir(Path(output_root).resolve())
    before = _snapshot_originals(project, output_root)
    staged = stage_baseline(project, output_root, waterfall_cache=waterfall_cache)
    private_json(staged / "original_files_before.json", before)
    log = staged / "baseline.log"
    module_name = "_abtrise_baseline_" + staged.name.replace("-", "_")
    old_path = list(sys.path)
    old_cwd = Path.cwd()
    old_umask = os.umask(0o077)
    os.environ.setdefault("MPLBACKEND", "Agg")
    try:
        # CLIFpy defaults also resolve against cwd, not only notebook __file__.
        os.chdir(staged)
        sys.path.insert(0, str(staged / "code"))
        spec = importlib.util.spec_from_file_location(module_name, staged / "code" / "01_cohort.py")
        if spec is None or spec.loader is None:
            raise RuntimeError("Cannot load the staged baseline notebook")
        notebook = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = notebook
        with log.open("w") as stream, redirect_stdout(stream), redirect_stderr(stream):
            spec.loader.exec_module(notebook)
            outputs, definitions = notebook.app.run()
            cohort = definitions["cohort_df"]
            replay, method_ties = replay_baseline(definitions["resp_wf"], definitions["hosp_final"], staged / "code" / "01_cohort.py", return_branch_candidates=True)
            parity = assert_baseline_parity(cohort, replay, method_ties=method_ties)
            private_parquet(staged / "baseline_method_ties.parquet", method_ties)
            panel = _eligible_panel(definitions)
            observations = (
                definitions["rs_all_df"]
                .join(panel.select("hospitalization_id").unique(), on="hospitalization_id", how="semi")
                .select([c for c in RESP_FIELDS if c in definitions["rs_all_df"].columns])
            )
            private_parquet(staged / "original_observations.parquet", observations)
            private_parquet(staged / "eligible_icu_windows.parquet", panel)
            private_parquet(staged / "baseline_events.parquet", replay)
            del observations, panel, definitions, outputs
            gc.collect()
    except Exception:
        # The detailed traceback stays in a private file, never in a public console.
        import traceback
        with log.open("a") as stream:
            traceback.print_exc(file=stream)
        raise RuntimeError(f"Staged baseline execution failed; inspect private log: {log}") from None
    finally:
        os.umask(old_umask)
        os.chdir(old_cwd)
        sys.path[:] = old_path
        sys.modules.pop(module_name, None)
        if log.exists() and os.name == "posix":
            log.chmod(0o600)
    after = _snapshot_originals(project, output_root)
    if before != after:
        private_json(staged / "original_files_after.json", after)
        raise RuntimeError("Original baseline files changed during staged execution")
    config = json.loads((staged / "clif_config.json").read_text())
    data_root = Path(config["data_directory"])
    respiratory_source = data_root / "clif_respiratory_support.parquet"
    metadata = {
        "staged_project": str(staged),
        "cohort": str(staged / "output_phi" / "cohort.parquet"),
        "waterfall": str(staged / "output_phi" / "resp_waterfall_cohort.parquet"),
        "original_observations": str(staged / "original_observations.parquet"),
        "eligible_windows": str(staged / "eligible_icu_windows.parquet"),
        "notebook_sha256": file_sha256(project / "code" / "01_cohort.py"),
        "source_parquet_sha256": file_sha256(respiratory_source),
        "original_files_unchanged": before == after,
        "parity": parity,
        "timezone_assumption": config.get("timezone"),
        "time_semantics": "Existing CLIF loader conversion followed by the unchanged notebook's strip_tz; source timezone interpretation not independently validated",
        "original_cache_copied": waterfall_cache is None and (project / "output_phi" / "resp_waterfall_cohort.parquet").is_file(),
        "explicit_cache_source": str(waterfall_cache) if waterfall_cache is not None else None,
        "waterfall_sha256": file_sha256(staged / "output_phi" / "resp_waterfall_cohort.parquet"),
        "versions": {name: importlib.metadata.version(name) for name in ("clifpy", "polars", "pandas", "marimo", "pyarrow")},
    }
    private_json(staged / "baseline_metadata.json", metadata)
    private_json(output_root / "latest_baseline.json", metadata)
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--waterfall-cache", type=Path, help="Explicit existing cache copied into staging; source is recorded")
    args = parser.parse_args()
    output_root = args.output_root or args.project / "output_phi" / "imv_comparison"
    try:
        result = run_baseline(args.project, output_root, waterfall_cache=args.waterfall_cache)
    except (ValueError, RuntimeError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 1
    print(json.dumps({"parity": result["parity"], "original_files_unchanged": result["original_files_unchanged"], "metadata": str(output_root / "latest_baseline.json")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
