"""ANCOM-BC-style bias-corrected CLR + Welch t (Lin & Peddada 2020, simplified).

This is a careful, documented simplification — not a reimplementation of the
R/Bioconductor package. See docs/methods/differential.md for scope notes.
"""
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def ancom_bc_da(
    counts: pd.DataFrame,
    labels: pd.Series,
    ref_group: str,
    alpha: float = 0.05,
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    labels = labels.loc[counts.index]
    logx = np.log(counts + pseudocount)
    # Per-taxon bias correction: subtract column median across samples
    bias = logx.median(axis=0)
    corrected = logx - bias

    ref_mask = (labels == ref_group).values
    alt_mask = ~ref_mask

    rows = []
    for taxon in corrected.columns:
        a = corrected.loc[ref_mask, taxon].values
        b = corrected.loc[alt_mask, taxon].values
        if a.std() == 0 and b.std() == 0:
            t, p = 0.0, 1.0
        else:
            t, p = ttest_ind(a, b, equal_var=False)
        rows.append(
            {
                "taxon": taxon,
                "t": float(t) if np.isfinite(t) else 0.0,
                "p": float(p) if np.isfinite(p) else 1.0,
                "effect": float(b.mean() - a.mean()),
            }
        )

    df = pd.DataFrame(rows)
    df["q"] = multipletests(df["p"].fillna(1.0), method="fdr_bh")[1]
    df["W_rank"] = df["effect"].abs() * -np.log10(df["q"].clip(lower=1e-300))
    df["significant"] = df["q"] < alpha
    return df.sort_values("W_rank", ascending=False).reset_index(drop=True)
