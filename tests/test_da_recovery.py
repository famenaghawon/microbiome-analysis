"""DA recovery + guard tests. The tests that prove the whole thing works."""
import numpy as np
from microbiome_analysis.data.synthetic import make_synthetic
from microbiome_analysis.preprocess.transforms import clr
from microbiome_analysis.differential.wilcoxon import wilcoxon_da
from microbiome_analysis.differential.ancom_bc import ancom_bc_da


def test_wilcoxon_recovers_planted_taxa():
    bundle = make_synthetic(seed=42, effect=3.0)
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
    bundle = make_synthetic(seed=42, effect=3.0)
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
    assert frac_sig < 0.15, f"Null FDR {frac_sig:.2%} too high"


def test_synthetic_actually_plants_effect():
    """Guard: if the generator silently fails to apply effect,
    planted taxa must have measurably higher mean abundance in group 1.

    NOTE: group labels are STRINGS ('0'/'1') - do not compare to int 0/1,
    that gives all-False mask, empty selection, and .mean() == nan.
    """
    b = make_synthetic(seed=42, effect=3.0)
    counts = b["counts"]
    group = b["metadata"]["group"]           # Series of str
    planted = b["truth"]["planted_taxa"]

    for taxon in planted:
        m0 = counts.loc[group == "0", taxon].mean()
        m1 = counts.loc[group == "1", taxon].mean()
        ratio = m1 / max(m0, 1e-9)
        assert ratio > 1.5, (
            f"Planting bug: {taxon} ratio {ratio:.2f} "
            f"(m0={m0:.1f}, m1={m1:.1f})"
        )


def test_synthetic_effect_absent_when_effect_is_one():
    """Complementary guard: with effect=1.0, no group difference exists."""
    b = make_synthetic(seed=42, effect=1.0)
    counts = b["counts"]
    group = b["metadata"]["group"]
    for taxon in b["truth"]["planted_taxa"]:
        m0 = counts.loc[group == "0", taxon].mean()
        m1 = counts.loc[group == "1", taxon].mean()
        ratio = m1 / max(m0, 1e-9)
        assert 0.7 < ratio < 1.4, f"{taxon} ratio {ratio:.2f} unexpected under null"
