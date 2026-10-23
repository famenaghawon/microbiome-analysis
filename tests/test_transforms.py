"""Compositional transforms + the depth bug guard."""
import numpy as np
import pandas as pd

from microbiome_analysis.preprocess.transforms import clr, relative, sample_depths


def test_depth_is_sum_not_mean_of_means():
    df = pd.DataFrame({"a": [10, 20], "b": [30, 40]}, index=["s1", "s2"])
    depths = sample_depths(df)
    assert depths["s1"] == 40
    assert depths["s2"] == 60
    # The classic bug would have produced 25.0 for both
    assert depths["s1"] != df.values.mean()


def test_relative_rows_sum_to_one():
    df = pd.DataFrame({"a": [10, 20], "b": [30, 40]}, index=["s1", "s2"])
    rel = relative(df)
    np.testing.assert_allclose(rel.sum(axis=1).values, [1.0, 1.0])


def test_clr_rows_sum_to_zero():
    df = pd.DataFrame({"a": [10, 20], "b": [30, 40]}, index=["s1", "s2"])
    c = clr(df)
    np.testing.assert_allclose(c.sum(axis=1).values, [0.0, 0.0], atol=1e-10)
