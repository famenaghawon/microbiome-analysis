"""Loader contract test. Uses a tiny hand-written fixture, no network."""
from pathlib import Path

import pandas as pd

from microbiome_analysis.data.loaders import load_hmp

FIX = Path("data/fixtures")


def _write_fixture():
    FIX.mkdir(parents=True, exist_ok=True)
    samples = [f"S{i:02d}" for i in range(6)]
    genera = [f"g__Taxon{i}" for i in range(5)]
    counts = pd.DataFrame(
        [[100, 50, 10, 5, 1]] * 6, index=samples, columns=genera
    )
    meta = pd.DataFrame(
        {"body_site": ["gut", "gut", "gut", "oral", "oral", "oral"]},
        index=samples,
    )
    counts.to_csv(FIX / "hmp_genus.tsv", sep="\t")
    meta.to_csv(FIX / "hmp_metadata.tsv", sep="\t")


def test_hmp_loader_reads_fixture():
    _write_fixture()
    bundle = load_hmp(path=str(FIX))
    assert bundle["source"] == "hmp"
    assert list(bundle["counts"].index) == [f"S{i:02d}" for i in range(6)]
    assert bundle["metadata"]["body_site"].tolist()[:3] == ["gut", "gut", "gut"]


def test_hmp_loader_rejects_low_join_rate():
    _write_fixture()
    # Corrupt metadata IDs to force mismatch
    meta_path = FIX / "hmp_metadata.tsv"
    df = pd.read_csv(meta_path, sep="\t", index_col=0)
    df.index = [f"X{i}" for i in range(len(df))]
    df.to_csv(meta_path, sep="\t")
    try:
        try:
            load_hmp(path=str(FIX))
            assert False, "should have raised"
        except ValueError as e:
            assert "Join rate" in str(e)
    finally:
        _write_fixture()  # restore
