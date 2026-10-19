"""The test that proves the whole thing works."""
import numpy as np
from microbiome_analysis.data.synthetic import make_synthetic
from microbiome_analysis.preprocess.transforms import clr
from microbiome_analysis.differential.wilcoxon import wilcoxon_da
from microbiome_analysis.differential.ancom_bc import ancom_bc_da


def test_wilcoxon_recovers_planted_taxa():
    bundle = make_synthetic(seed=42)
    result = wilcoxon_da(
        clr(bundle["counts"]),
        bundle["metadata"]["group"],
        group_a="0",
        group_b="1",
    )
    top_k = len(bundle["truth"]["planted_taxa"])
    recovered = set(result.head(top_k)["taxon"]) & set(bundle["truth"]["planted_taxa"])
    recall = len(recovered) / top_k
    assert recall >= 0.6, f"Wilcoxon recall {recall:.2f} < 0.6"
    assert (result["q"] < 0.05).sum() >= top_k * 0.5


def test_ancom_bc_recovers_planted_taxa():
    bundle = make_synthetic(seed=42)
    result = ancom_bc_da(
        bundle["counts"],
        bundle["metadata"]["group"],
        ref_group="0",
    )
    top_k = len(bundle["truth"]["planted_taxa"])
    recovered = set(result.head(top_k)["taxon"]) & set(bundle["truth"]["planted_taxa"])
    recall = len(recovered) / top_k
    assert recall >= 0.6, f"ANCOM-BC recall {recall:.2f} < 0.6"


def test_fdr_controls_false_positives():
    """Under the null (effect=1.0), few taxa should be called significant."""
    bundle = make_synthetic(seed=7, effect=1.0)
    result = wilcoxon_da(
        clr(bundle["counts"]),
        bundle["metadata"]["group"],
        group_a="0",
        group_b="1",
    )
    frac_sig = (result["q"] < 0.05).mean()
    assert frac_sig < 0.15, f"Null FDR {frac_sig:.2%} too high — check FDR handling"
