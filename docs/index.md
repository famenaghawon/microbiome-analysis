# microbiome-analysis

> Reproducible microbiome analysis on real HMP-adjacent data:
> compositionality-aware differential abundance, diversity, ordination,
> and leakage-safe ML. Audited by tests and CI.

## What this project does

End-to-end pipeline that runs from **real HMP-adjacent 16S data**
(122 samples x 122 genera) and produces:

- **Differential abundance** with BH-FDR (Wilcoxon on CLR + ANCOM-BC-style)
- **Diversity analysis** (richness, Shannon, Simpson; Bray-Curtis + Jaccard)
- **Ordination** (PCoA on beta distance matrices)
- **PERMANOVA** (group-difference test on community composition)
- **Leakage-safe ML** (balanced accuracy, AUC, permutation-ranked features)

## Quick start

```powershell
pip install -e ".[dev,dashboard]"
pytest -q
python scripts/run_pipeline.py --source gmrepo --out results/real_sanity/
python scripts/run_ml.py
streamlit run dashboard/app.py
