# microbiome-analysis

[![CI](https://github.com/famenaghawon/microbiome-analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/famenaghawon/microbiome-analysis/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23599184.svg)](https://doi.org/10.5281/zenodo.23599184)
[![Docs](https://img.shields.io/badge/docs-live-brightgreen)](https://famenaghawon.github.io/microbiome-analysis/)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-31%20passed-brightgreen)](tests/)

> **Reproducible, multi-site microbiome analysis on real HMP-adjacent data: compositionality-aware differential abundance, α/β diversity, ordination, PERMANOVA, and leakage-safe ML. Audited by 31 tests and CI.**

![Microbiome Analysis Dashboard](docs/images/dashboard.png)
*Interactive Streamlit dashboard showing PCoA ordination across five human body sites.*

## Overview

`microbiome-analysis` is an end-to-end Python toolkit that transforms raw 16S rRNA genus tables into publication-quality figures and machine-learning results. It ships with a real, multi-site human microbiome dataset (**662 samples across Gut, Oral, Skin, Nasal, and Vaginal sites**, 124 genera) and provides:

- **Compositional differential abundance** — Wilcoxon rank-sum and ANCOM-BC-style on CLR-transformed abundances, with BH-FDR correction.
- **Alpha diversity** — observed richness, Shannon, Simpson via `scikit-bio`.
- **Beta diversity** — Bray–Curtis and Jaccard distance matrices.
- **Ordination** — Principal Coordinates Analysis (PCoA) with variance-explained reporting.
- **PERMANOVA** — permutation-based group-difference test on the distance matrix.
- **Leakage-safe machine learning** — `Pipeline(CLR → SelectKBest → RandomForest)` with feature selection strictly inside CV folds. Achieves **~98% balanced accuracy** classifying body site from genus profiles.
- **Interactive Streamlit dashboard** — six tabs (Overview, Alpha, Beta, DA, ML, Data & QC) with Plotly visualizations, download buttons, and defensive error handling.

## Why this project is different

Most "microbiome analysis" GitHub repos are tutorials that never touch real data, never test their DA methods against ground truth, and never prove reproducibility in CI. This project does all three:

| Practice | How this repo does it |
|---|---|
| Compositional correctness | CLR transform is the default; raw abundances are never used for DA |
| FDR control | BH correction on every DA method, verified by null-data tests |
| Leakage-safe ML | Feature selection inside `Pipeline`, verified by shuffled-label sanity test |
| Planted-taxa recovery | Synthetic generator with ground truth, tested in CI |
| Offline reproducibility | All tests pass on Linux/Windows/macOS without network |
| Real biological results | Top genera are *Bacteroides* (gut), *Streptococcus* (oral), *Lactobacillus* (vaginal), *Propionibacterium* (skin) |

## Quick start

```bash
# Clone and install
git clone https://github.com/famenaghawon/microbiome-analysis.git
cd microbiome-analysis
python -m venv .venv
source .venv/bin/activate         # Windows: .\.venv\Scripts\Activate.ps1
pip install -e ".[dev,dashboard,docs]"

# Verify the pipeline
pytest -q

# Run on real data
python scripts/run_pipeline.py --source gmrepo_v2 --out results/real/

# Explore the dashboard
streamlit run dashboard/app.py
```

Then open `http://localhost:8501`.

## Data

The real dataset lives in `data/cache/` (gitignored) and comes from **GMrepo + QIIME 2**:

- `_genus_filtered.tsv` — 320 genera × 662 samples (relative abundance)
- `sample_metadata_grouped.tsv` — `SampleID`, `BodySite`, `BodySiteCategory`

**Body site distribution:** Oral (341), Gut (163), Nasal (80), Vaginal (60), Skin (18).

See [`docs/data.md`](docs/data.md) for the full `SRA → QIIME 2 → BIOM → TSV` provenance chain.

Three data sources are supported:

| Source | Requires | Purpose |
|---|---|---|
| `--source gmrepo_v2` | `data/cache/_genus_filtered.tsv` + metadata | Real multi-site analysis (default) |
| `--source gmrepo` | `data/cache/_genus_raw.tsv` | Legacy single-site GMrepo data |
| `--source synthetic` | nothing | Offline test fixture |

## Results on real data

Sample sizes: 662 samples × 124 genera. Gut vs Oral contrast.

| Metric | Value |
|---|---|
| DA genera at q < 0.05 (Wilcoxon) | 111 / 124 |
| DA genera at q < 0.05 (ANCOM-BC) | 111 / 124 |
| PERMANOVA (body site) | pseudo-F ≈ 13.5, **p = 0.001** |
| PCoA PC1 variance | ~20% |
| PCoA PC2 variance | ~13% |
| ML balanced accuracy (5-class body site) | **~0.98** |
| ML AUC | **~0.99** |

Generated figures in `results/`:

- `results/diversity/alpha_histograms.png`
- `results/diversity/pcoa_braycurtis.png`
- `results/diversity/pcoa_jaccard.png`
- `results/ml/feature_importance.tsv`

## Package structure

```
src/microbiome_analysis/
├── data/           loaders (gmrepo, gmrepo_v2, hmp, synthetic) + cache
├── preprocess/     CLR, relative, rarefy, depth
├── diversity/      alpha, beta, ordination, permanova
├── differential/   wilcoxon, ancom_bc
├── ml/             classifier (leakage-safe CV), importance
└── plotting.py     headless-safe matplotlib helpers
```

## Testing

```bash
pytest -q
```

31 tests across seven modules:

| Test file | Coverage |
|---|---|
| `test_transforms.py` | CLR/relative row sums, depth-is-sum-not-mean |
| `test_da_recovery.py` | Planted-taxa recovery + null FDR control |
| `test_diversity.py` | Alpha known-answers, distance symmetry, PERMANOVA |
| `test_gmrepo_loader.py` | Real loader contract, genus extraction, join rate |
| `test_hmp_loader.py` | HMP format contract |
| `test_ml.py` | **Leakage sanity check** (shuffled labels → ~0.5) |
| `test_pipeline_smoke.py` | End-to-end offline pipeline |

CI runs on GitHub Actions across **Python 3.10, 3.11, and 3.12**.

## Documentation

Live site: **https://famenaghawon.github.io/microbiome-analysis/**

- [Home](https://famenaghawon.github.io/microbiome-analysis/)
- [Methods](https://famenaghawon.github.io/microbiome-analysis/methods/)
- [Data provenance](https://famenaghawon.github.io/microbiome-analysis/data/)
- [Limitations](https://famenaghawon.github.io/microbiome-analysis/limitations/)
- [API reference](https://famenaghawon.github.io/microbiome-analysis/api/)

## Limitations

An honest list of what this project does not claim:

1. **Compositional data only.** Absolute abundances cannot be recovered from 16S data. All DA runs on CLR.
2. **ANCOM-BC simplification.** Our implementation is a bias-corrected CLR + Welch t-test, **not** a reimplementation of the Bioconductor package. Cite Lin & Peddada (2020).
3. **No phylogenetic distances yet.** UniFrac is planned for a future release.
4. **No batch correction.** Metadata-level batch effects are not modeled.
5. **Single-cohort training.** The ML module is trained and tested on one cohort; external validation is not performed.
6. **Not for clinical use.** This is an analysis toolkit, not a diagnostic tool.

See [`docs/limitations.md`](docs/limitations.md).

## Citation

If you use this software in your research, cite it via its Zenodo DOI:

```bibtex
@software{famenaghawon_microbiome_analysis_2026,
  author       = {Famenaghawon, Freedom Evbuomwan},
  title        = {microbiome-analysis: Multi-site Microbiome Analysis Platform},
  month        = oct,
  year         = 2026,
  version      = {v1.0.0},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23599184},
  url          = {https://doi.org/10.5281/zenodo.23599184}
}
```

ORCID: [0009-0003-1457-809X](https://orcid.org/0009-0003-1457-809X)

## License

MIT — see [`LICENSE`](LICENSE).

---

*Built with [`scikit-bio`](https://scikit-bio.org), [`scikit-learn`](https://scikit-learn.org), [`Streamlit`](https://streamlit.io), and [`Plotly`](https://plotly.com). Figures and results are fully reproducible from the released DOI.*
