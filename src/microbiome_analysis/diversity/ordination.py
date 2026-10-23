"""PCoA ordination for beta diversity distances.

Principal Coordinates Analysis (PCoA) = MDS on a distance matrix.
Reduces N x N distances to 2D coordinates for plotting.
"""
import numpy as np
import pandas as pd
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa


def ordination_coords(
    distances: DistanceMatrix, n_components: int = 2
) -> tuple[pd.DataFrame, np.ndarray]:
    """Run PCoA, return sample coordinates + explained variance ratio.

    Args:
        distances: N x N DistanceMatrix
        n_components: number of coordinates to keep (default 2 for plotting)

    Returns:
        (coords_df, explained_ratio)
        coords_df: DataFrame indexed by sample_id, columns PC1, PC2, ...
        explained_ratio: array of variance explained by each PC (sums to 1)
    """
    result = pcoa(distances)
    coords = result.samples.iloc[:, :n_components].copy()
    coords.index.name = "sample_id"
    coords.columns = [f"PC{i + 1}" for i in range(n_components)]
    explained = np.array(result.proportion_explained[:n_components])
    return coords, explained
