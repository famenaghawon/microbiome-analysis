"""One-command pipeline entry point.

Sources:
    --source gmrepo      Real GMrepo/QIIME2 genus abundance table (default)
    --source hmp        Real HMP-format cache (hmp_genus.tsv + hmp_metadata.tsv)
    --source synthetic  Seeded synthetic with planted DA taxa (tests only)

Examples:
    python scripts/run_pipeline.py --source gmrepo --out results/real/
    python scripts/run_pipeline.py --source gmrepo --group-col depth_bin \
        --group-a low --group-b high --out results/real/
    python scripts/run_pipeline.py --source synthetic --seed 42 --out results/synthetic/

Produces under --out:
    summary.json       Full provenance + QC + top taxa + DA summary
    wilcoxon_da.tsv    Per-taxon Wilcoxon results with BH-FDR
    ancom_bc_da.tsv    Per-taxon ANCOM-BC-style results
    counts_top50.tsv   Top 50 genera by total reads (for manual inspection)
    depths.tsv         Per-sample sequencing depth
"""
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from microbiome_analysis.data.loaders import load
from microbiome_analysis.data.synthetic import make_synthetic
from microbiome_analysis.differential.ancom_bc import ancom_bc_da
from microbiome_analysis.differential.wilcoxon import wilcoxon_da
from microbiome_analysis.preprocess.transforms import clr, sample_depths


def _group_labels(meta, group_col, ga, gb):
    """Return (labels_series, group_a, group_b, contrast_mode)."""
    if group_col == "depth_bin":
        if "depth" not in meta.columns:
            raise ValueError("depth_bin needs meta['depth'] -- internal bug")
        med = meta["depth"].median()
        labels = pd.Series(
            ["low" if d < med else "high" for d in meta["depth"]], index=meta.index
        )
        return labels, "low", "high", "depth_bin"

    if group_col not in meta.columns:
        raise SystemExit(
            f"'{group_col}' not in metadata. Available columns: {list(meta.columns)}"
        )
    if ga is None or gb is None:
        groups = sorted(meta[group_col].dropna().unique().tolist())
        if len(groups) < 2:
            raise SystemExit(
                f"'{group_col}' has only {groups}; need 2 groups for DA. "
                "Pass --group-col depth_bin for a synthetic sanity check."
            )
        ga = ga or groups[0]
        gb = gb or groups[1]
    return meta[group_col], ga, gb, "metadata"


def main():
    ap = argparse.ArgumentParser(
        description="Microbiome analysis pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "--source",
        default="gmrepo",
        choices=["gmrepo", "hmp", "synthetic"],
        help="Data source (default: gmrepo -- real data)",
    )
    ap.add_argument("--path", default=None, help="Override path for the source")
    ap.add_argument("--max-samples", type=int, default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--group-col",
        default="depth_bin",
        help="Metadata column for DA contrast. 'depth_bin' (default) = sanity "
        "split; use e.g. 'disease' or 'body_site' for real biology.",
    )
    ap.add_argument("--group-a", default=None)
    ap.add_argument("--group-b", default=None)
    ap.add_argument("--out", default="results")
    ap.add_argument(
        "--transform",
        default="clr",
        choices=["clr", "relative"],
        help="Abundance transform before Wilcoxon DA (default: clr)",
    )
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"[pipeline] loading source={args.source}", file=sys.stderr)
    if args.source == "synthetic":
        bundle = make_synthetic(seed=args.seed)
        if args.group_col == "depth_bin":
            args.group_col = "group"
    else:
        bundle = load(args.source, path=args.path,
                      max_samples=args.max_samples, seed=args.seed)

    counts = bundle["counts"]
    meta = bundle["metadata"].copy()

    if counts.shape[0] < 4:
        raise SystemExit(
            f"Only {counts.shape[0]} samples after QC -- too few for DA."
        )

    depths = sample_depths(counts)
    meta["depth"] = depths.values

    if args.transform == "clr":
        x = clr(counts)
    else:
        x = counts.div(counts.sum(axis=1), axis=0)

    labels, ga, gb, mode = _group_labels(
        meta, args.group_col, args.group_a, args.group_b
    )
    print(
        f"[pipeline] DA contrast: {ga!r} vs {gb!r} "
        f"(mode={mode}, n_a={(labels == ga).sum()}, n_b={(labels == gb).sum()})",
        file=sys.stderr,
    )

    print("[pipeline] running Wilcoxon DA...", file=sys.stderr)
    w = wilcoxon_da(x, labels, group_a=ga, group_b=gb)

    print("[pipeline] running ANCOM-BC-style DA...", file=sys.stderr)
    a = ancom_bc_da(counts, labels, ref_group=ga)

    totals = counts.sum().sort_values(ascending=False)
    top50 = totals.head(50).to_frame("total_reads")
    top50["rel_abundance"] = top50["total_reads"] / depths.sum()

    bundle_report = bundle.get("report", {})
    summary = {
        "run": {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "source": args.source,
            "path": args.path or bundle_report.get("path", "n/a"),
            "seed": args.seed,
            "transform": args.transform,
        },
        "data": {
            "n_samples": int(counts.shape[0]),
            "n_genera": int(counts.shape[1]),
            "median_depth": float(depths.median()),
            "min_depth": int(depths.min()),
            "max_depth": int(depths.max()),
            "total_reads": int(depths.sum()),
        },
        "loader_report": bundle_report,
        "contrast": {
            "group_col": args.group_col,
            "group_a": ga,
            "group_b": gb,
            "n_a": int((labels == ga).sum()),
            "n_b": int((labels == gb).sum()),
            "mode": mode,
            "note": (
                "depth_bin = sanity split, NOT biology"
                if mode == "depth_bin"
                else "user-specified metadata contrast"
            ),
        },
        "top_genera": (
            top50.head(20)
            .reset_index()
            .rename(columns=lambda c: "taxon" if c in ("index", "__genus__", "taxonomy", "sample_id") else c)
            .to_dict(orient="records")
        ),
        "wilcoxon": {
            "n_significant_q05": int(w["significant"].sum()),
            "top20": w.head(20)[["taxon", "effect", "q", "prev_a", "prev_b"]].to_dict(
                orient="records"
            ),
        },
        "ancom_bc": {
            "n_significant_q05": int(a["significant"].sum()),
            "top20": a.head(20)[["taxon", "effect", "q"]].to_dict(orient="records"),
        },
    }

    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    w.to_csv(out / "wilcoxon_da.tsv", sep="\t", index=False)
    a.to_csv(out / "ancom_bc_da.tsv", sep="\t", index=False)
    top50.to_csv(out / "counts_top50.tsv", sep="\t")
    depths.to_frame().to_csv(out / "depths.tsv", sep="\t")

    print(
        json.dumps(
            {
                "source": args.source,
                "n_samples": summary["data"]["n_samples"],
                "n_genera": summary["data"]["n_genera"],
                "median_depth": summary["data"]["median_depth"],
                "contrast": summary["contrast"],
                "wilcoxon_sig_q05": summary["wilcoxon"]["n_significant_q05"],
                "top_genera_total_reads": [
                    r["taxon"] for r in summary["top_genera"][:5]
                ],
            },
            indent=2,
        )
    )
    print(f"[pipeline] wrote results to {out.resolve()}", file=sys.stderr)


if __name__ == "__main__":
    main()
