"""Leakage-safe ML classifier on microbiome profiles.

Critical design principle: feature selection MUST happen inside CV folds,
never on the full dataset before CV. Otherwise the model "peeks" at test
data, inflating reported accuracy. This module enforces that via sklearn
Pipeline, which refits every step on each train fold.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer

from microbiome_analysis.preprocess.transforms import clr


def build_pipeline(k: int = 50, seed: int = 0) -> Pipeline:
    """Return sklearn Pipeline: CLR -> SelectKBest -> RandomForest.

    Fitting this pipeline on train-fold-only data guarantees no leakage:
    the selector sees only train-fold features when picking top-k.
    """
    return Pipeline(
        [
            ("clr", FunctionTransformer(clr)),
            ("select", SelectKBest(f_classif, k=k)),
            ("clf", RandomForestClassifier(
                n_estimators=500,
                max_depth=None,
                min_samples_leaf=2,
                random_state=seed,
                n_jobs=-1,
            )),
        ]
    )


def evaluate(
    counts: pd.DataFrame,
    labels: pd.Series,
    k: int = 50,
    folds: int = 5,
    seed: int = 0,
) -> dict:
    """Cross-validated evaluation, leakage-safe.

    Returns dict with:
        balanced_accuracy_mean, balanced_accuracy_std
        auc_mean, auc_std
        per_fold_balanced_accuracy
        per_fold_auc
        confusion_matrix_last_fold
    """
    labels = labels.loc[counts.index]
    y_encoded, classes = pd.factorize(labels)
    n_classes = len(classes)
    if n_classes < 2:
        raise ValueError(
            f"Need >= 2 classes for classification. Got: {list(classes)}"
        )

    cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    ba_scores, auc_scores = [], []
    last_cm = None

    for train_idx, test_idx in cv.split(counts, y_encoded):
        X_train = counts.iloc[train_idx]
        X_test = counts.iloc[test_idx]
        y_train = y_encoded[train_idx]
        y_test = y_encoded[test_idx]

        pipe = build_pipeline(k=k, seed=seed)
        pipe.fit(X_train, y_train)
        y_pred = pipe.predict(X_test)

        ba_scores.append(balanced_accuracy_score(y_test, y_pred))
        last_cm = confusion_matrix(y_test, y_pred, labels=list(range(n_classes)))

        if n_classes == 2:
            proba = pipe.predict_proba(X_test)[:, 1]
            try:
                auc_scores.append(roc_auc_score(y_test, proba))
            except ValueError:
                pass

    result = {
        "balanced_accuracy_mean": float(np.mean(ba_scores)),
        "balanced_accuracy_std": float(np.std(ba_scores)),
        "per_fold_balanced_accuracy": [float(x) for x in ba_scores],
        "n_classes": n_classes,
        "classes": [str(c) for c in classes],
        "confusion_matrix_last_fold": last_cm.tolist() if last_cm is not None else None,
        "k_features": k,
        "folds": folds,
    }
    if auc_scores:
        result["auc_mean"] = float(np.mean(auc_scores))
        result["auc_std"] = float(np.std(auc_scores))
        result["per_fold_auc"] = [float(x) for x in auc_scores]
    return result


def train_final_and_rank_features(
    counts: pd.DataFrame,
    labels: pd.Series,
    k: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """Fit on all data (for final model) and rank top-k features by importance.

    Returns DataFrame with taxon + importance, sorted desc.
    Useful for interpretability: "which genera matter most for prediction?"
    """
    labels = labels.loc[counts.index]
    y_encoded, _classes = pd.factorize(labels)

    pipe = build_pipeline(k=k, seed=seed)
    pipe.fit(counts, y_encoded)

    selector = pipe.named_steps["select"]
    clf = pipe.named_steps["clf"]

    mask = selector.get_support()
    selected_taxa = counts.columns[mask]
    importances = clf.feature_importances_

    df = pd.DataFrame(
        {"taxon": selected_taxa, "importance": importances}
    ).sort_values("importance", ascending=False)
    df["rank"] = range(1, len(df) + 1)
    return df.reset_index(drop=True)


def leakage_sanity_check(seed: int = 7, k: int = 30) -> dict:
    """If the pipeline leaks, this returns balanced_accuracy >> 0.5 on shuffled labels.

    For binary classification with 5-fold CV, chance is 0.5.
    A leakage-free pipeline trained on shuffled labels must be ~0.5.
    """
    from microbiome_analysis.data.synthetic import make_synthetic

    bundle = make_synthetic(seed=seed)
    counts = bundle["counts"]
    rng = np.random.default_rng(0)
    fake_labels = pd.Series(
        rng.integers(0, 2, len(counts)).astype(str), index=counts.index
    )
    result = evaluate(counts, fake_labels, k=k, seed=seed)
    return result
