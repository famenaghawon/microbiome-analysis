"""Alpha diversity metrics via scikit-bio.

Alpha diversity = "how diverse is each sample individually?"
Common metrics:
    - richness (observed_otus): count of distinct taxa
    - Shannon: entropy-based, weighs both richness and evenness
    - Simpson: probability two random reads are different taxa
    - Faith's PD: phylogenetic diversity (needs a tree, not implemented yet)
"""
import pandas as pd
from skbio.diversity import alpha_diversity


def alpha_table(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample alpha diversity.

    Args:
        counts_df: samples x taxa, integer read counts.

    Returns:
        DataFrame indexed by sample_id, with columns:
        richness, shannon, simpson.
    """
    values = counts_df.values
    ids = [str(s) for s in counts_df.index]

    richness = alpha_diversity("observed_otus", values, ids=ids)
    shannon = alpha_diversity("shannon", values, ids=ids)
    simpson = alpha_diversity("simpson", values, ids=ids)

    return pd.DataFrame(
        {
            "richness": richness.values if hasattr(richness, "values") else richness,
            "shannon": shannon.values if hasattr(shannon, "values") else shannon,
            "simpson": simpson.values if hasattr(simpson, "values") else simpson,
        },
        index=pd.Index(ids, name="sample_id"),
    )
