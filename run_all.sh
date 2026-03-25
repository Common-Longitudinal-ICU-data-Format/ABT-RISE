#!/usr/bin/env bash
set -e

uv run python code/01_cohort.py
uv run python code/02_wide_dataset.py
uv run python code/03_sat.py
uv run python code/04_sbt_both.py
uv run python code/05_analysis_dataset.py
