# ml/classifier.py
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier

def build_pipeline(k=50, seed=0):
    return Pipeline([
        ("clr",           FunctionTransformer(clr)),        # fit-free transform
        ("select",        SelectKBest(f_classif, k=k)),      # fit on TRAIN folds only
        ("clf",           RandomForestClassifier(500, random_state=seed)),
    ])

def evaluate(counts, labels, folds=5, seed=0):
    cv = StratifiedKFold(folds, shuffle=True, random_state=seed)
    res = cross_validate(build_pipeline(seed=seed), counts, labels, cv=cv,
                         scoring=["balanced_accuracy", "roc_auc_ovr_weighted"])
    return {k: (v.mean(), v.std()) for k, v in res.items() if k.startswith("test_")}

def test_no_feature_selection_on_test():
    """Sanity: if a leakage-free pipeline is fit on shuffled labels,
    balanced accuracy should be ~chance."""
    counts, truth = make_synthetic(seed=7)
    rng = np.random.default_rng(0)
    fake_labels = rng.integers(0, 2, len(counts))
    m = evaluate(counts, fake_labels)
    assert m["test_balanced_accuracy"][0] < 0.65
