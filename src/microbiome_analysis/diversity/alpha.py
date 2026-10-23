"""Alpha diversity metrics via scikit-bio."""
import pandas as pd
from skbio.diversity import alpha_diversity


def alpha_table(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample richness, Shannon, and Simpson diversity."""
    values = counts_df.values
    ids = list(counts_df.index)
    return pd.DataFrame(
        {
            "richness": alpha_diversity("observed_otus", values, ids=ids),
            "shannon": alpha_diversity("shannon", values, ids=ids),
            "simpson": alpha_diversity("simpson", values, ids=ids),
        }
    )
