"""Loader tests for real GMrepo data. Fixture is a 20x10 subset."""
from pathlib import Path
import pandas as pd
import pytest
from microbiome_analysis.data.gmrepo import load_gmrepo, _extract_genus


FIXTURE = Path("data/fixtures/gmrepo_mini.tsv")


def test_extract_genus_strips_full_taxonomy():
    s = (
        "d__Bacteria;k__Bacillati;p__Bacteroidota;c__Bacteroidia;"
        "o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
    )
    assert _extract_genus(s) == "g__Bacteroides"


def test_extract_genus_strips_other_suffix():
    s = "d__Bacteria;...;g__Actinomyces--other"
    assert _extract_genus(s) == "g__Actinomyces"


def test_extract_genus_rejects_uninformative():
    assert _extract_genus("d__Bacteria;p__Firmicutes;__;__;__;__;__") is None
    assert _extract_genus("d__Bacteria;g__uncultured") is None


def test_load_gmrepo_reads_fixture():
    if not FIXTURE.exists():
        pytest.skip("Fixture missing: " + str(FIXTURE))
    bundle = load_gmrepo(path=FIXTURE, min_depth=1, min_prevalence=0.0)
    assert bundle["source"] == "gmrepo"
    assert bundle["counts"].shape[1] > 0
    assert (bundle["counts"].dtypes == int).all()
    assert bundle["metadata"]["body_site"].iloc[0] == "gut"


def test_load_gmrepo_drops_low_depth_samples():
    if not FIXTURE.exists():
        pytest.skip("Fixture missing")
    raw = pd.read_csv(FIXTURE, sep="\t", index_col=0)
    dead = raw.columns[0]
    raw[dead] = 0
    tmp = Path("data/fixtures/_test_dead.tsv")
    raw.to_csv(tmp, sep="\t")
    try:
        bundle = load_gmrepo(path=tmp, min_depth=1, min_prevalence=0.0)
        assert dead not in bundle["counts"].index
    finally:
        tmp.unlink()
