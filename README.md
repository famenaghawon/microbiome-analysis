# microbiome-analysis

[![CI](https://github.com/USER/microbiome-analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/USER/microbiome-analysis/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Reproducible microbiome analysis from public data to diversity, differential abundance, and interpretable ML — audited by tests and CI.**

Runs end-to-end on a cached public subset (`--source hmp`) or a seeded synthetic cohort
with planted DA taxa (`--source synthetic`). Compositional-aware (CLR), FDR-controlled,
leakage-safe ML, one-command pipeline, Streamlit dashboard, and a CI job that re-runs
the whole thing on every PR.

## Honest scope

This is a reproducible analysis toolkit and educational reference — **not a clinical tool**.
See [docs/limitations.md](docs/limitations.md).

## Quick start

```bash
pip install -e ".[dev,dashboard]"
pytest -q
python scripts/run_pipeline.py --source synthetic --seed 42 --out results/
