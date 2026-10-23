"""PERMANOVA test for group differences in beta diversity.

PERMANOVA = multivariate ANOVA on distances.
Tests whether samples in different metadata groups have different
community composition, using permutation-based p-values.
"""
from skbio.stats.distance import DistanceMatrix, permanova


def run_permanova(
    distances: DistanceMatrix, groups, permutations: int = 999
) -> dict:
    """Run PERMANOVA on a distance matrix grouped by metadata.

    Args:
        distances: N x N DistanceMatrix with sample IDs
        groups: pd.Series indexed by sample ID, values = group labels
        permutations: number of permutations for p-value

    Returns:
        dict with method, pseudo_f, p_value, n_samples, n_groups
    """
    # reindex aligns groups to distance-matrix sample order.
    # NOTE: do NOT use .loc[list] here -- pandas treats list as column
    # axis and raises IndexingError on a Series. .reindex is the correct API.
    aligned = groups.reindex(distances.ids)

    # skbio.permanova returns a pandas Series with string keys
    result = permanova(distances, aligned.values, permutations=permutations)

    return {
        "method": "PERMANOVA",
        "pseudo_f": float(result["test statistic"]),
        "p_value": float(result["p-value"]),
        "n_samples": int(result["sample size"]),
        "n_groups": int(result["number of groups"]),
        "permutations": permutations,
    }
