"""Beta diversity distance matrices via scikit-bio.

Beta diversity = "how different are two samples from each other?"
Produces an N x N distance matrix (N = number of samples).
Common metrics:
    - Bray-Curtis: abundance-aware, most common for 16S
    - Jaccard: presence/absence only
    - UniFrac: phylogeny-aware (needs a tree, not implemented yet)
"""
import pandas as pd
from skbio.stats.distance import DistanceMatrix


def beta_distances(counts_df: pd.DataFrame) -> dict:
    """Compute Bray-Curtis and Jaccard distance matrices.

    Args:
        counts_df: samples x taxa, integer read counts.

    Returns:
        dict of {"braycurtis": DistanceMatrix, "jaccard": DistanceMatrix}
        Each DistanceMatrix is N x N with sample IDs as labels.
    """
    from skbio.diversity import beta_diversity

    values = counts_df.values
    ids = [str(s) for s in counts_df.index]

    bray = beta_diversity("braycurtis", values, ids=ids)
    jaccard = beta_diversity("jaccard", values, ids=ids)

    return {
        "braycurtis": bray,
        "jaccard": jaccard,
    }


def distances_to_frame(distances: DistanceMatrix) -> pd.DataFrame:
    """Convert a DistanceMatrix to a pandas DataFrame (for easy inspection)."""
    return pd.DataFrame(distances.data, index=distances.ids, columns=distances.ids)
