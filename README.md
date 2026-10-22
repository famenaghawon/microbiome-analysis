# microbiome-analysis

[![CI](https://github.com/USER/microbiome-analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/USER/microbiome-analysis/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Reproducible microbiome analysis on real HMP-adjacent data: compositionality-aware
> differential abundance, FDR control, planted-taxa recovery tests, and a one-command
> pipeline that ships from raw 16S to results.**

## What this is

A working, tested Python pipeline that runs end-to-end on **real HMP-adjacent
16S data** (122 samples, 122 genera, real gut microbiome from QIIME 2 /
DADA2 / Greengenes2). It computes:

- Compositional transforms (CLR, relative abundance, rarefaction)
- Differential abundance (Wilcoxon on CLR, ANCOM-BC-style, with BH-FDR)
- **Verified pipeline correctness** via synthetic planted-taxa recovery
- **Real biological signal** — top genera are Bacteroides, Faecalibacterium,
  Alistipes, Agathobacter — the canonical healthy gut microbiome

## Honest scope

A reproducible analysis toolkit and educational reference -- **not** a
clinical tool. See [`docs/limitations.md`](docs/limitations.md).

## Quick start

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[dev,dashboard]"
pytest -q
python scripts/run_pipeline.py --source gmrepo --out results/real_sanity/
