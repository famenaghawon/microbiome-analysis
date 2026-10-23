"""Loader for the new cleaned GMrepo dataset (genus × sample TSV + metadata TSV).

Format:
    _genus_filtered.tsv:
        - Rows: genus names (already cleaned, e.g. "Bacteroides")
        - Columns: sample IDs (strings, e.g. "700014386")
        - Values: relative abundances (sum to 1.0 per sample)

    sample_metadata_grouped.tsv:
        - Column: SampleID (matches abundance columns)
        - Column: BodySite (e.g. "Stool", "Oral", "Skin")
        - Column: BodySiteCategory (e.g. "Gut", "Oral", "Skin")

Returns the standard bundle:
    {
        "counts":   DataFrame (samples × genera),
        "metadata": DataFrame (samples × [BodySite, BodySiteCategory]),
        "source":   "gmrepo_v2",
        "report":   dict of QC info
    }
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ABUNDANCE_PATH = Path("data/cache/_genus_filtered.tsv")
META_PATH = Path("data/cache/sample_metadata_grouped.tsv")


def load_gmrepo_v2(
    abundance_path: str | Path | None = None,
    meta_path: str | Path | None = None,
    min_prevalence: float = 0.05,
    min_rel_abundance: float = 1e-6,
):
    ab_path = Path(abundance_path) if abundance_path else ABUNDANCE_PATH
    mt_path = Path(meta_path) if meta_path else META_PATH

    if not ab_path.exists():
        raise FileNotFoundError(f"Abundance table not found: {ab_path}")
    if not mt_path.exists():
        raise FileNotFoundError(f"Metadata not found: {mt_path}")

    print(f"[gmrepo_v2] loading {ab_path.name}", file=sys.stderr)
    # Read: rows = genera, columns = samples
    raw = pd.read_csv(ab_path, sep="\t", index_col=0)
    n_genus, n_sample = raw.shape
    print(f"[gmrepo_v2] raw: {n_genus} genera × {n_sample} samples", file=sys.stderr)

    # Transpose: samples × genera
    counts = raw.T
    counts.index.name = "sample_id"
    # CRITICAL: force string dtype to match metadata (avoid int vs str join miss)
    counts.index = counts.index.astype(str)

    # Sanity: check if values look like relative abundances (sum to ~1 per row)
    row_sums = counts.sum(axis=1)
    median_sum = float(row_sums.median())
    print(
        f"[gmrepo_v2] median row sum: {median_sum:.4f} "
        f"({'relative abundance' if median_sum < 10 else 'counts'})",
        file=sys.stderr,
    )
    is_relative = median_sum < 10

    # Read metadata
    print(f"[gmrepo_v2] loading {mt_path.name}", file=sys.stderr)
    meta = pd.read_csv(mt_path, sep="\t")
    # Find the sample ID column (could be SampleID, sample_id, etc.)
    id_col = None
    for candidate in ["SampleID", "sample_id", "Sample", "sample"]:
        if candidate in meta.columns:
            id_col = candidate
            break
    if id_col is None:
        raise ValueError(
            f"No sample ID column found. Columns: {list(meta.columns)}"
        )
    meta = meta.set_index(id_col)
    meta.index.name = "sample_id"
    meta.index = meta.index.astype(str) 
    print(
        f"[gmrepo_v2] metadata: {len(meta)} rows, "
        f"columns: {list(meta.columns)}",
        file=sys.stderr,
    )

    # Join: keep only samples in both
    inter = counts.index.intersection(meta.index)
    print(
        f"[gmrepo_v2] joined: {len(inter)} samples "
        f"(dropped {len(counts) - len(inter)} table-only, "
        f"{len(meta) - len(inter)} meta-only)",
        file=sys.stderr,
    )
    if len(inter) < len(counts) * 0.9:
        raise ValueError(
            f"Join rate {len(inter) / len(counts):.2%} < 90% — sample IDs don't match"
        )
    counts = counts.loc[inter]
    meta = meta.loc[inter]

    # Filter: drop rare genera (low prevalence)
    prevalence = (counts > 0).mean(axis=0)
    keep = prevalence >= min_prevalence
    n_dropped = int((~keep).sum())
    counts = counts.loc[:, keep]
    print(
        f"[gmrepo_v2] dropped {n_dropped} genera with "
        f"prevalence < {min_prevalence:.0%}",
        file=sys.stderr,
    )

    # Filter: drop very low-abundance genera
    mean_abund = counts.mean(axis=0)
    keep2 = mean_abund >= min_rel_abundance
    n_dropped2 = int((~keep2).sum())
    counts = counts.loc[:, keep2]
    print(
        f"[gmrepo_v2] dropped {n_dropped2} genera with "
        f"mean abundance < {min_rel_abundance}",
        file=sys.stderr,
    )

    report = {
        "source": "gmrepo_v2",
        "abundance_path": str(ab_path),
        "metadata_path": str(mt_path),
        "n_raw_genera": int(n_genus),
        "n_raw_samples": int(n_sample),
        "n_dropped_low_prevalence": int(n_dropped),
        "n_dropped_low_abundance": int(n_dropped2),
        "n_final_samples": int(counts.shape[0]),
        "n_final_genera": int(counts.shape[1]),
        "median_row_sum": float(row_sums.median()),
        "is_relative_abundance": bool(is_relative),
        "body_site_categories": (
            dict(meta["BodySiteCategory"].value_counts())
            if "BodySiteCategory" in meta.columns else {}
        ),
    }

    print(
        f"[gmrepo_v2] final: {counts.shape[0]} samples × {counts.shape[1]} genera",
        file=sys.stderr,
    )
    print(
        f"[gmrepo_v2] body site categories: {report['body_site_categories']}",
        file=sys.stderr,
    )

    return {
        "counts": counts,
        "metadata": meta,
        "source": "gmrepo_v2",
        "report": report,
    }
