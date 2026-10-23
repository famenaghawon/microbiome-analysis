"""Leakage-safe ML classifier tests."""
import pandas as pd

from microbiome_analysis.data.synthetic import make_synthetic
from microbiome_analysis.ml.classifier import (
    build_pipeline,
    evaluate,
    leakage_sanity_check,
    train_final_and_rank_features,
)


def test_leakage_sanity_chance_accuracy():
    """Shuffled labels must yield balanced accuracy near 0.5 (chance)."""
    result = leakage_sanity_check(seed=7, k=30)
    ba = result["balanced_accuracy_mean"]
    assert 0.35 < ba < 0.65, (
        f"Leakage suspected: balanced accuracy {ba} on shuffled labels. "
        "Expected ~0.5 for leakage-free binary classification."
    )


def test_pipeline_classifies_planted_signal():
    """Planted DA signal must be learnable above chance."""
    bundle = make_synthetic(seed=42, effect=3.0)
    counts = bundle["counts"]
    labels = pd.Series(bundle["metadata"]["group"].values, index=counts.index)
    result = evaluate(counts, labels, k=30, folds=5, seed=0)
    assert result["balanced_accuracy_mean"] > 0.7, (
        "Expected > 0.7 on planted signal, got {}".format(
            result["balanced_accuracy_mean"]
        )
    )
    assert result["n_classes"] == 2


def test_rank_features_returns_valid_table():
    bundle = make_synthetic(seed=42, effect=3.0)
    counts = bundle["counts"]
    labels = pd.Series(bundle["metadata"]["group"].values, index=counts.index)
    df = train_final_and_rank_features(counts, labels, k=20, seed=0)
    assert len(df) == 20
    assert set(df.columns) == {"taxon", "importance", "rank"}
    # Importance is non-negative and sorted
    assert (df["importance"] >= 0).all()
    assert df["importance"].is_monotonic_decreasing
    # Rank starts at 1
    assert df["rank"].iloc[0] == 1


def test_pipeline_builds_without_error():
    pipe = build_pipeline(k=10, seed=42)
    assert hasattr(pipe, "fit")
    assert hasattr(pipe, "predict")
