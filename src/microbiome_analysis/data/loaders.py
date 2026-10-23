"""Data loaders for real and synthetic microbiome data.

Resolution for --source gmrepo:
    data/cache/_genus_raw.tsv      (real QIIME2/GMrepo wide-format table)

Resolution for --source hmp:
    1. Explicit path (--path) -- must exist
    2. data/cache/hmp_genus.tsv + data/cache/hmp_metadata.tsv
    3. data/fixtures/hmp_genus.tsv + data/fixtures/hmp_metadata.tsv
    4. $env:MB_HMP_URL (base URL hosting both files)
    5. Refuse loudly -- never silently fall back to synthetic.

For --source synthetic:
    Seeded generator with planted DA taxa (test-only).

All loaders return:
    {
        "counts":   DataFrame (samples x genera), integer counts
        "metadata": DataFrame (samples x ...), at least 'body_site'
        "source":   str
        "report":   dict of QC info (optional)
    }
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

from .cache import CACHE_DIR, fetch
from .synthetic import make_synthetic

HMP_TABLE_NAME = "hmp_genus.tsv"
HMP_META_NAME = "hmp_metadata.tsv"
FIXTURES_DIR = Path("data/fixtures")
GMREPO_DEFAULT = Path("data/cache/_genus_raw.tsv")


def load_gmrepo(
    path=None,
    min_depth=1000,
    min_prevalence=0.10,
    max_samples=None,
    seed=42,
):
    """Load real QIIME2/GMrepo wide-format genus abundance table."""
    from .gmrepo import load_gmrepo as _load_gmrepo

    p = Path(path) if path else GMREPO_DEFAULT
    bundle = _load_gmrepo(
        path=p, min_depth=min_depth, min_prevalence=min_prevalence
    )

    if max_samples and bundle["counts"].shape[0] > max_samples:
        keep = bundle["counts"].sample(max_samples, random_state=seed).index
        bundle["counts"] = bundle["counts"].loc[keep]
        bundle["metadata"] = bundle["metadata"].loc[keep]
        bundle["report"]["n_final_samples"] = len(keep)
        bundle["report"]["max_samples_applied"] = int(max_samples)
        print(
            f"[gmrepo] subsampled to {max_samples} samples (seed={seed})",
            file=sys.stderr,
        )

    return bundle


def load_gmrepo_dispatch(**kwargs):
    """Convenience wrapper for dispatcher tests."""
    return load_gmrepo(**kwargs)


def _validate_hmp(counts, meta):
    """Join on sample ID, filter to valid body sites, print join diagnostics."""
    if "body_site" not in meta.columns:
        raise ValueError("metadata missing required column 'body_site'")
    meta = meta.copy()
    meta["body_site"] = meta["body_site"].astype(str).str.lower().str.strip()

    inter = counts.index.intersection(meta.index)
    drop_table = len(counts) - len(inter)
    drop_meta = len(meta) - len(inter)
    join_rate = len(inter) / len(counts) if len(counts) else 0.0

    print(
        f"[hmp] joined {len(inter)} samples (dropped {drop_table} table-only, {drop_meta} meta-only; "
        f"join rate {join_rate:.2%})",
        file=sys.stderr,
    )
    if join_rate < 0.8:
        raise ValueError(
            f"Join rate {join_rate:.2%} < 80% -- sample IDs don't line up. "
            f"First table IDs: {list(counts.index[:3])}, first meta IDs: {list(meta.index[:3])}"
        )

    counts = counts.loc[inter]
    meta = meta.loc[inter]
    valid_sites = {"gut", "oral", "skin", "nasal", "vaginal",
                   "stool", "saliva", "tongue"}
    mask = meta["body_site"].isin(valid_sites)
    if mask.sum() == 0:
        raise ValueError(
            "None of the body_site values are recognized. "
            "Found: {}. Expected one of: {}".format(
                sorted(meta["body_site"].unique())[:10], sorted(valid_sites))
        )
    meta = meta[mask]
    counts = counts.loc[meta.index]
    print(
        "[hmp] after body-site filter: {} samples ({})".format(
            len(counts), dict(meta["body_site"].value_counts())),
        file=sys.stderr,
    )
    return counts, meta


def load_hmp(
    path=None,
    max_samples=None,
    seed=42,
    allow_synthetic_fallback=False,
):
    """Load real HMP-format genus table + metadata. Refuses silent synthetic fallback."""
    counts_path = meta_path = None

    if path:
        p = Path(path)
        if p.is_file() and p.name == HMP_TABLE_NAME:
            counts_path, meta_path = p, p.parent / HMP_META_NAME
        elif p.is_dir():
            counts_path, meta_path = p / HMP_TABLE_NAME, p / HMP_META_NAME
        if not (counts_path and counts_path.exists()
                and meta_path and meta_path.exists()):
            raise FileNotFoundError(
                f"Expected {HMP_TABLE_NAME} and {HMP_META_NAME} under {path}"
            )
    elif (CACHE_DIR / HMP_TABLE_NAME).exists() and (CACHE_DIR / HMP_META_NAME).exists():
        counts_path = CACHE_DIR / HMP_TABLE_NAME
        meta_path = CACHE_DIR / HMP_META_NAME
    elif (FIXTURES_DIR / HMP_TABLE_NAME).exists() and (FIXTURES_DIR / HMP_META_NAME).exists():
        counts_path = FIXTURES_DIR / HMP_TABLE_NAME
        meta_path = FIXTURES_DIR / HMP_META_NAME
    elif os.getenv("MB_HMP_URL"):
        base = os.environ["MB_HMP_URL"].rstrip("/")
        counts_path = fetch(f"{base}/{HMP_TABLE_NAME}",
                            CACHE_DIR / HMP_TABLE_NAME)
        meta_path = fetch(f"{base}/{HMP_META_NAME}",
                          CACHE_DIR / HMP_META_NAME)
    else:
        msg = (
            "No real HMP data available.\n"
            f"  Put {HMP_TABLE_NAME} and {HMP_META_NAME} in {CACHE_DIR.resolve()}/, OR\n"
            f"  Put them in {FIXTURES_DIR.resolve()}/, OR\n"
            "  Set $env:MB_HMP_URL to a base URL hosting both files."
        )
        if allow_synthetic_fallback:
            print("[hmp] WARNING: falling back to synthetic (tests only).\n" + msg,
                  file=sys.stderr)
            return make_synthetic(seed=seed)
        raise FileNotFoundError(msg)

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    counts, meta = _validate_hmp(counts, meta)

    if max_samples and len(counts) > max_samples:
        counts = counts.sample(max_samples, random_state=seed)
        meta = meta.loc[counts.index]

    return {
        "counts": counts,
        "metadata": meta,
        "source": "hmp",
        "report": {
            "source": "hmp",
            "counts_path": str(counts_path),
            "meta_path": str(meta_path),
            "n_samples": len(counts),
            "n_genera": int(counts.shape[1]),
        },
    }


def load(source, path=None, max_samples=None, seed=42):
    """Dispatch to the correct loader by source name."""
    if source == "gmrepo_v2":
        from .gmrepo_v2 import load_gmrepo_v2
        bundle = load_gmrepo_v2()
        if max_samples and bundle["counts"].shape[0] > max_samples:
            keep = bundle["counts"].sample(max_samples, random_state=seed).index
            bundle["counts"] = bundle["counts"].loc[keep]
            bundle["metadata"] = bundle["metadata"].loc[keep]
            bundle["report"]["n_final_samples"] = int(len(keep))
            bundle["report"]["max_samples_applied"] = int(max_samples)
        return bundle
    if source == "gmrepo":
        return load_gmrepo(path=path, max_samples=max_samples, seed=seed)
    if source == "hmp":
        return load_hmp(path=path, max_samples=max_samples, seed=seed)
    if source == "synthetic":
        return make_synthetic(seed=seed)
    raise ValueError(
        "Unknown source: {!r}. Expected one of: 'gmrepo_v2', 'gmrepo', 'hmp', 'synthetic'".format(source)
    )
