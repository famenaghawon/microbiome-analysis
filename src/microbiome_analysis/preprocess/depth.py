# src/microbiome_analysis/preprocess/depth.py
def sample_depths(counts_df) -> pd.Series:
    """Per-sample sequencing depth. NEVER mean().mean()."""
    return counts_df.sum(axis=1).rename("depth")

