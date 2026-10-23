"""End-to-end smoke test. Runs the full pipeline offline on a fixture.

This test does NOT require network access, does NOT require the real
cache file, and does NOT require a specific working directory. It uses
a 20x10 subset of the real HMP-adjacent data stored at
data/fixtures/gmrepo_mini.tsv.
"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

FIXTURE = Path("data/fixtures/gmrepo_mini.tsv")


@pytest.fixture(scope="module")
def ensure_fixture():
    """Extract a small fixture from the real cache if it's not there yet."""
    if FIXTURE.exists():
        return FIXTURE

    real = Path("data/cache/_genus_raw.tsv")
    if not real.exists():
        pytest.skip("Real cache not present -- smoke test needs fixture or cache")

    raw = pd.read_csv(real, sep="\t", index_col=0)
    mini = raw.iloc[:20, :10]
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    mini.to_csv(FIXTURE, sep="\t")
    return FIXTURE


def test_loaders_import():
    """Both loaders must import cleanly."""
    from microbiome_analysis.data.gmrepo import load_gmrepo  # noqa: F401
    from microbiome_analysis.data.loaders import (  # noqa: F401
        load,
        load_hmp,
    )


def test_synthetic_loader_contract():
    """make_synthetic must return the standard bundle."""
    from microbiome_analysis.data.synthetic import make_synthetic

    b = make_synthetic(seed=42)
    assert b["source"] == "synthetic"
    assert b["counts"].shape == (100, 80)
    assert "report" in b
    assert "truth" in b
    assert len(b["truth"]["planted_taxa"]) == 10


def test_gmrepo_loader_on_fixture(ensure_fixture):
    """Loader reads the fixture, extracts genera, aggregates correctly."""
    from microbiome_analysis.data.gmrepo import load_gmrepo

    b = load_gmrepo(path=ensure_fixture, min_depth=1, min_prevalence=0.0)
    assert b["source"] == "gmrepo"
    assert b["counts"].shape[1] > 0
    assert all(str(c).startswith("g__") for c in b["counts"].columns)
    assert (b["counts"].dtypes == int).all()
    assert b["metadata"]["body_site"].iloc[0] == "gut"


def test_pipeline_runs_offline_on_fixture(ensure_fixture, tmp_path):
    """The whole pipeline runs end-to-end and produces expected files."""
    out = tmp_path / "smoke_out"
    result = subprocess.run(
        [
            sys.executable,
            "scripts/run_pipeline.py",
            "--source", "gmrepo",
            "--path", str(ensure_fixture),
            "--out", str(out),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        "Pipeline failed:\nSTDOUT:\n" + result.stdout +
        "\nSTDERR:\n" + result.stderr
    )

    for fname in ("summary.json", "wilcoxon_da.tsv",
                  "ancom_bc_da.tsv", "counts_top50.tsv", "depths.tsv"):
        p = out / fname
        assert p.exists(), "Missing output: " + fname
        assert p.stat().st_size > 0, "Empty output: " + fname

    summary = json.loads((out / "summary.json").read_text())
    assert summary["run"]["source"] == "gmrepo"
    assert summary["data"]["n_samples"] > 0
    assert summary["data"]["n_genera"] > 0
    assert "top_genera" in summary
    assert "wilcoxon" in summary
    assert "ancom_bc" in summary
    assert "contrast" in summary


def test_pipeline_synthetic_offline(tmp_path):
    """The pipeline also works on synthetic data (regression guard)."""
    out = tmp_path / "smoke_synth"
    result = subprocess.run(
        [
            sys.executable,
            "scripts/run_pipeline.py",
            "--source", "synthetic",
            "--seed", "42",
            "--out", str(out),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    summary = json.loads((out / "summary.json").read_text())
    assert summary["run"]["source"] == "synthetic"
    assert summary["data"]["n_samples"] == 100
