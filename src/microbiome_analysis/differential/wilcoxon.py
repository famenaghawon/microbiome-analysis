"""Wilcoxon rank-sum DA with BH-FDR — the honest baseline."""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def wilcoxon_da(
    clr_df: pd.DataFrame,
    labels: pd.Series,
    group_a: str,
    group_b: str,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Per-taxon Welch-style Wilcoxon on CLR abundances, BH-FDR corrected.

    Returns columns: taxon, statistic, p, q, effect (mean_b - mean_a),
    prev_a, prev_b, direction, significant.
    """
    labels = labels.loc[clr_df.index]
    a_mask = (labels == group_a).values
    b_mask = (labels == group_b).values
    if a_mask.sum() < 3 or b_mask.sum() < 3:
        raise ValueError(f"Need >=3 samples per group. Got a={a_mask.sum()}, b={b_mask.sum()}")

    rows = []
    for taxon in clr_df.columns:
        a = clr_df.loc[a_mask, taxon].values
        b = clr_df.loc[b_mask, taxon].values
        try:
            stat, p = mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            stat, p = np.nan, 1.0
        rows.append(
            {
                "taxon": taxon,
                "statistic": float(stat) if not np.isnan(stat) else 0.0,
                "p": float(p),
                "effect": float(np.mean(b) - np.mean(a)),
                "prev_a": float((a > 0).mean()),
                "prev_b": float((b > 0).mean()),
            }
        )

    df = pd.DataFrame(rows)
    df["q"] = multipletests(df["p"].fillna(1.0), method="fdr_bh")[1]
    df["direction"] = np.where(df["effect"] > 0, f"up_in_{group_b}", f"up_in_{group_a}")
    df["significant"] = df["q"] < alpha
    return df.sort_values(["q", "effect"], ascending=[True, False]).reset_index(drop=True)
