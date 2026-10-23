"""Per-sample sequencing depth. NEVER mean().mean()."""
import pandas as pd


def sample_depths(counts_df) -> pd.Series:
    """Per-sample sequencing depth. NEVER mean().mean()."""
    return counts_df.sum(axis=1).rename("depth")
