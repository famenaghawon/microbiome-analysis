"""Loader for real GMrepo / QIIME2 genus abundance tables.

Input format (tab-separated, wide):
    taxonomy                         SRR41063392  SRR41063393  ...
    d__Bacteria;k__...;g__Bacteroides    1561.0        2048.0    ...

- Column 0: full taxonomic string
- Columns 1..N: sample IDs (one per column)
- Values: integer read counts

Output:
    {
        "counts":   DataFrame (samples x genera), integer counts
        "metadata": DataFrame (samples x ...), body_site / group placeholder
        "source":   "gmrepo"
        "report":   dict of QC info
    }
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

GENUS_RE = re.compile(r";g__([^;]+)")
EXCLUDED_GENUS = {"__", "__other", "Incertae_Sedis", "uncultured", "other"}


def _extract_genus(taxonomy: str):
    """'d__Bacteria;...;g__Bacteroides' -> 'g__Bacteroides'. None if uninformative."""
    m = GENUS_RE.search(taxonomy)
    if not m:
        return None
    name = m.group(1).strip()
    if not name or name in EXCLUDED_GENUS or name.startswith("__"):
        return None
    name = name.split("--")[0]
    return f"g__{name}"


def load_gmrepo(
    path="data/cache/_genus_raw.tsv",
    min_depth=1000,
    min_prevalence=0.10,
):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Real GMrepo table not found at {path}.\n"
            "Place the wide-format genus table at data/cache/_genus_raw.tsv.\n"
            "Format: col 0 = taxonomy string, cols 1..N = samples, values = counts."
            
        )

    raw = pd.read_csv(path, sep="\t", index_col=0)
    n_raw_taxa, n_raw_samples = raw.shape
    print(f"[gmrepo] loaded {n_raw_taxa} taxa x {n_raw_samples} samples",
          file=sys.stderr)

    genus_map = pd.Series([_extract_genus(t) for t in raw.index], index=raw.index)
    keep_mask = genus_map.notna()
    n_dropped_taxa = int((~keep_mask).sum())

    raw = raw.loc[keep_mask]
    genus_map = genus_map.loc[keep_mask]

    raw["__genus__"] = genus_map.values
    aggregated = raw.groupby("__genus__", sort=False).sum()
    n_unique_genera = aggregated.shape[0]

    counts = aggregated.T.astype(int)
    counts.index.name = "sample_id"

    depths = counts.sum(axis=1)
    dead = depths < min_depth
    n_dead = int(dead.sum())
    if n_dead:
        print(f"[gmrepo] dropping {n_dead} sample(s) with depth < {min_depth}: {depths[dead].to_dict()}", file=sys.stderr)
    counts = counts.loc[~dead]

    prevalence = (counts > 0).mean(axis=0)
    rare = prevalence < min_prevalence
    n_rare = int(rare.sum())
    if n_rare:
        print(f"[gmrepo] dropping {n_rare} taxa with prevalence < {min_prevalence:.0%}", file=sys.stderr)
    counts = counts.loc[:, ~rare]

    metadata = pd.DataFrame(
        {"body_site": "gut", "group": "all"}, index=counts.index
    )

    report = {
        "source": "gmrepo",
        "path": str(path),
        "n_raw_taxa": int(n_raw_taxa),
        "n_raw_samples": int(n_raw_samples),
        "n_unique_genera": int(n_unique_genera),
        "n_dropped_taxa_uninformative": int(n_dropped_taxa),
        "n_dropped_taxa_rare": int(n_rare),
        "n_dropped_samples_dead": int(n_dead),
        "n_final_samples": int(counts.shape[0]),
        "n_final_genera": int(counts.shape[1]),
        "median_depth": float(depths[~dead].median()),
        "min_depth": int(depths[~dead].min()),
        "max_depth": int(depths[~dead].max()),
    }
    print("[gmrepo] final: {} samples x {} genera (median depth {:.0f})".format(
        counts.shape[0], counts.shape[1], report["median_depth"]), file=sys.stderr)

    return {"counts": counts, "metadata": metadata,
            "source": "gmrepo", "report": report}
