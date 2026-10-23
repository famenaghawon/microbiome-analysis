"""Leakage-safe ML classifier on microbiome profiles."""
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer

from microbiome_analysis.preprocess.transforms import clr


def build_pipeline(k: int = 50, seed: int = 0) -> Pipeline:
    """Return a scikit-learn pipeline that transforms then classifies."""
    return Pipeline(
        [
            ("clr", FunctionTransformer(clr)),
            ("select", SelectKBest(f_classif, k=k)),
            ("clf", RandomForestClassifier(500, random_state=seed)),
        ]
    )


def evaluate(counts, labels, folds: int = 5, seed: int = 0):
    """Cross-validated evaluation with leakage-safe feature selection."""
    cv = StratifiedKFold(folds, shuffle=True, random_state=seed)
    res = cross_validate(
        build_pipeline(seed=seed),
        counts,
        labels,
        cv=cv,
        scoring=["balanced_accuracy", "roc_auc_ovr_weighted"],
    )
    return {
        k: (v.mean(), v.std()) for k, v in res.items() if k.startswith("test_")
    }


def sanity_check_leakage():
    """If the pipeline leaks, shuffled labels should give ~chance accuracy."""
    from microbiome_analysis.data.synthetic import make_synthetic

    bundle = make_synthetic(seed=7)
    counts = bundle["counts"]
    rng = np.random.default_rng(0)
    fake_labels = rng.integers(0, 2, len(counts))
    return evaluate(counts, fake_labels)
