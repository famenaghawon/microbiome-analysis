"""Diversity, ordination, PERMANOVA tests with known-answer fixtures."""
import numpy as np
import pandas as pd
import pytest

from microbiome_analysis.diversity.alpha import alpha_table
from microbiome_analysis.diversity.beta import beta_distances
from microbiome_analysis.diversity.ordination import ordination_coords
from microbiome_analysis.diversity.permanova import run_permanova


@pytest.fixture
def simple_counts():
    """Two samples, one uniform, one dominated. Known-answer for tests."""
    return pd.DataFrame(
        {
            "t1": [10, 90],
            "t2": [10, 5],
            "t3": [10, 5],
            "t4": [10, 0],
        },
        index=["uniform", "dominated"],
    )


def test_alpha_uniform_is_highest_shannon(simple_counts):
    alpha = alpha_table(simple_counts)
    assert alpha.loc["uniform", "shannon"] > alpha.loc["dominated", "shannon"]


def test_alpha_richness_equals_nonzero_count(simple_counts):
    alpha = alpha_table(simple_counts)
    assert alpha.loc["uniform", "richness"] == 4
    assert alpha.loc["dominated", "richness"] == 3


def test_alpha_schema(simple_counts):
    alpha = alpha_table(simple_counts)
    assert list(alpha.index) == ["uniform", "dominated"]
    assert set(alpha.columns) == {"richness", "shannon", "simpson"}


def test_beta_symmetric_zero_diagonal(simple_counts):
    dists = beta_distances(simple_counts)
    bray = dists["braycurtis"].data
    # Symmetric
    assert np.allclose(bray, bray.T)
    # Zero diagonal
    assert np.allclose(np.diag(bray), 0)


def test_beta_range(simple_counts):
    dists = beta_distances(simple_counts)
    bray = dists["braycurtis"].data
    assert bray.min() >= 0
    assert bray.max() <= 1


def test_pcoa_output_shape(simple_counts):
    dists = beta_distances(simple_counts)["braycurtis"]
    coords, explained = ordination_coords(dists)
    assert coords.shape[0] == 2  # two samples
    assert coords.shape[1] == 2  # two components
    assert len(explained) == 2


def test_permanova_runs_and_returns_schema():
    """PERMANOVA returns the expected dict schema on a small fixture."""
    # 4 samples x 4 taxa, two clearly distinct community types
    counts = pd.DataFrame(
        {
            "t1": [100, 0, 90, 5],
            "t2": [0, 100, 5, 90],
            "t3": [50, 50, 10, 10],
            "t4": [0, 0, 0, 0],
        },
        index=["A1", "A2", "B1", "B2"],
    )
    dists = beta_distances(counts)["braycurtis"]
    groups = pd.Series(["g1", "g1", "g2", "g2"], index=["A1", "A2", "B1", "B2"])
    result = run_permanova(dists, groups, permutations=99)
    assert result["method"] == "PERMANOVA"
    assert 0 <= result["p_value"] <= 1
    assert result["n_samples"] == 4
    assert result["n_groups"] == 2
    assert result["permutations"] == 99

