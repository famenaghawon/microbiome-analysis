"""Real HMP genus-level loader with honest provenance.

Resolution order for --source hmp:
  1. Explicit path (--path) — wins, must exist.
  2. Local cache file: data/cache/hmp_genus.tsv (+ hmp_metadata.tsv).
  3. $MB_HMP_URL (you host the real table on Zenodo/HF and set this).
  4. Refuse loudly — never silently fall back to synthetic.

Format contract for hmp_genus.tsv (tab-separated):
  - First column = sample ID (index).
  - Remaining columns = genus names, prefixed `g__` (e.g. g__Bacteroides).
  - Values = integer read counts.
Format contract for hmp_metadata.tsv:
  - First column = sample ID (must match table index).
  - Column `body_site` ∈ {gut, oral, skin, nasal, vaginal} (case-insensitive;
    the loader normalizes to lowercase).
"""
import os
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

from .cache import CACHE_DIR, fetch
from .synthetic import make_synthetic


HMP_TABLE_NAME = "hmp_genus.tsv"
HMP_META_NAME = "hmp_metadata.tsv"


def _validate(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Join, drop non-intersecting samples, log what happened."""
    if "body_site" not in meta.columns:
        raise ValueError("metadata missing required column 'body_site'")
    meta = meta.copy()
    meta["body_site"] = meta["body_site"].astype(str).str.lower().str.strip()

    inter = counts.index.intersection(meta.index)
    drop_table = len(counts) - len(inter)
    drop_meta = len(meta) - len(inter)
    join_rate = len(inter) / len(counts) if len(counts) else 0.0
    print(
        f"[hmp] joined {len(inter)} samples "
        f"(dropped {drop_table} table-only, {drop_meta} meta-only; "
        f"join rate {join_rate:.2%})",
        file=sys.stderr,
    )
    if join_rate < 0.8:
        raise ValueError(
            f"Join rate {join_rate:.2%} < 80% — sample IDs don't line up. "
            f"Check first IDs: table={list(counts.index[:3])}, "
            f"meta={list(meta.index[:3])}"
        )
    counts = counts.loc[inter]
    meta = meta.loc[inter]
    meta = meta[meta["body_site"].isin({"gut", "oral", "skin", "nasal", "vaginal"})]
    counts = counts.loc[meta.index]
    print(f"[hmp] after body-site filter: {len(counts)} samples", file=sys.stderr)
    return counts, meta


def load_hmp(
    path: Optional[str] = None,
    max_samples: Optional[int] = None,
    seed: int = 42,
    allow_synthetic_fallback: bool = False,
) -> dict:
    """Load real HMP genus table + metadata. Refuses silent synthetic fallback."""
    counts_path = meta_path = None

    if path:
        p = Path(path)
        if p.is_file() and p.name == HMP_TABLE_NAME:
            counts_path, meta_path = p, p.parent / HMP_META_NAME
        elif p.is_dir():
            counts_path, meta_path = p / HMP_TABLE_NAME, p / HMP_META_NAME
        if not (counts_path and counts_path.exists() and meta_path.exists()):
            raise FileNotFoundError(
                f"Expected {HMP_TABLE_NAME} and {HMP_META_NAME} under {path}"
            )
    elif (CACHE_DIR / HMP_TABLE_NAME).exists() and (CACHE_DIR / HMP_META_NAME).exists():
        counts_path = CACHE_DIR / HMP_TABLE_NAME
        meta_path = CACHE_DIR / HMP_META_NAME
    elif os.getenv("MB_HMP_URL"):
        base = os.environ["MB_HMP_URL"].rstrip("/")
        counts_path = fetch(f"{base}/{HMP_TABLE_NAME}", CACHE_DIR / HMP_TABLE_NAME)
        meta_path = fetch(f"{base}/{HMP_META_NAME}", CACHE_DIR / HMP_META_NAME)
    else:
        msg = (
            "No real HMP data available.\n"
            f"  Put {HMP_TABLE_NAME} and {HMP_META_NAME} in {CACHE_DIR.resolve()}/, OR\n"
            "  set $env:MB_HMP_URL to a base URL hosting both files.\n"
            "See docs/data.md for the format contract and suggested sources."
        )
        if allow_synthetic_fallback:
            print("[hmp] WARNING: falling back to synthetic (tests only).\n" + msg,
                  file=sys.stderr)
            return make_synthetic(seed=seed)
        raise FileNotFoundError(msg)

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    counts, meta = _validate(counts, meta)

    if max_samples and len(counts) > max_samples:
        counts = counts.sample(max_samples, random_state=seed)
        meta = meta.loc[counts.index]

    return {"counts": counts, "metadata": meta, "source": "hmp"}


def load(source: str, path: Optional[str] = None, max_samples: Optional[int] = None,
         seed: int = 42) -> dict:
    if source == "hmp":
        return load_hmp(path=path, max_samples=max_samples, seed=seed)
    if source == "synthetic":
        return make_synthetic(seed=seed)
    raise ValueError(f"Unknown source: {source}")
