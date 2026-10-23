"""Train leakage-safe classifier on real data; save metrics + feature importance."""
import json
from pathlib import Path

import pandas as pd

from microbiome_analysis.data.loaders import load
from microbiome_analysis.ml.classifier import evaluate, train_final_and_rank_features


def main():
    out = Path("results/ml")
    out.mkdir(parents=True, exist_ok=True)

    bundle = load("gmrepo")
    counts = bundle["counts"]
    print(f"Loaded {counts.shape[0]} samples x {counts.shape[1]} genera")

    # Target: depth_bin (median split of sequencing depth).
    # This is a genuine prediction task: can the model learn abundance
    # structure that correlates with sequencing depth?
    depths = counts.sum(axis=1)
    med = depths.median()
    labels = pd.Series(
        ["high" if d >= med else "low" for d in depths],
        index=counts.index,
    )
    print("Target distribution:", dict(labels.value_counts()))

    # Leakage-safe cross-validated evaluation
    print("Running 5-fold stratified CV...")
    metrics = evaluate(counts, labels, k=50, folds=5, seed=42)

    # Final model on all data, rank top 20 features
    print("Fitting final model, ranking top features...")
    importance = train_final_and_rank_features(counts, labels, k=20, seed=42)
    importance.to_csv(out / "feature_importance.tsv", sep="\t", index=False)

    summary = {
        "n_samples": int(counts.shape[0]),
        "n_genera": int(counts.shape[1]),
        "target": "depth_bin (median split)",
        "classes": metrics["classes"],
        "metrics": metrics,
        "top_10_features": importance.head(10)[["rank", "taxon", "importance"]].to_dict(
            orient="records"
        ),
    }
    (out / "metrics.json").write_text(json.dumps(summary, indent=2))
    print()
    print(json.dumps({k: v for k, v in summary.items() if k != "metrics"}, indent=2))
    print()
    print("CV balanced accuracy: {:.3f} +/- {:.3f}".format(
        metrics["balanced_accuracy_mean"], metrics["balanced_accuracy_std"]))
    if "auc_mean" in metrics:
        print("CV AUC: {:.3f} +/- {:.3f}".format(
            metrics["auc_mean"], metrics["auc_std"]))
    print(f"Wrote: {out.resolve()}")


if __name__ == "__main__":
    main()
