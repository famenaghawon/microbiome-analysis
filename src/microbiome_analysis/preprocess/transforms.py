"""Compositional transforms. Depth is NEVER computed via mean().mean()."""
import numpy as np
import pandas as pd


def relative(counts: pd.DataFrame) -> pd.DataFrame:
    return counts.div(counts.sum(axis=1), axis=0)


def clr(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    """Centered log-ratio transform. Samples x taxa in, same shape out."""
    x = np.log(counts + pseudocount)
    return x.sub(x.mean(axis=1), axis=0)


def rarefy(counts: pd.DataFrame, depth: int, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    out = np.zeros_like(counts.values, dtype=float)
    for i, row in enumerate(counts.values):
        total = row.sum()
        if total <= 0:
            continue
        p = row / total
        n = min(depth, int(total))
        out[i] = rng.multinomial(n, p)
    return pd.DataFrame(out, index=counts.index, columns=counts.columns)


def sample_depths(counts: pd.DataFrame) -> pd.Series:
    """Per-sample depth. The metric that MUST be sum, not mean-of-means."""
    return counts.sum(axis=1).rename("depth")
