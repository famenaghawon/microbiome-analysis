"""One-command pipeline entry point."""
import argparse
import json
from pathlib import Path

from microbiome_analysis.data.synthetic import make_synthetic
from microbiome_analysis.preprocess.transforms import clr, sample_depths
from microbiome_analysis.differential.wilcoxon import wilcoxon_da
from microbiome_analysis.differential.ancom_bc import ancom_bc_da


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="synthetic", choices=["synthetic", "hmp", "local"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.source != "synthetic":
        raise NotImplementedError(f"source={args.source} is Phase 2 — use 'synthetic' for now")

    bundle = make_synthetic(seed=args.seed)
    counts, meta = bundle["counts"], bundle["metadata"]
    depths = sample_depths(counts)

    c = clr(counts)
    w = wilcoxon_da(c, meta["group"], group_a="0", group_b="1")
    a = ancom_bc_da(counts, meta["group"], ref_group="0")

    summary = {
        "source": args.source,
        "seed": args.seed,
        "n_samples": int(len(counts)),
        "n_taxa": int(counts.shape[1]),
        "median_depth": float(depths.median()),
        "min_depth": int(depths.min()),
        "truth_planted": bundle["truth"]["planted_taxa"],
        "wilcoxon": {
            "n_significant": int(w["significant"].sum()),
            "top10": w.head(10)[["taxon", "effect", "q"]].to_dict(orient="records"),
        },
        "ancom_bc": {
            "n_significant": int(a["significant"].sum()),
            "top10": a.head(10)[["taxon", "effect", "q"]].to_dict(orient="records"),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    w.to_csv(out / "wilcoxon_da.tsv", sep="\t", index=False)
    a.to_csv(out / "ancom_bc_da.tsv", sep="\t", index=False)
    print(json.dumps({k: summary[k] for k in ("n_samples", "n_taxa", "median_depth")}, indent=2))
    print(f"Wrote results to {out.resolve()}")


if __name__ == "__main__":
    main()
