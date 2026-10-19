"""Synthetic generator with planted DA taxa. Ground truth for CI tests."""
import numpy as np
import pandas as pd


def make_synthetic(
    n_samples: int = 100,
    n_taxa: int = 80,
    planted: int = 10,
    effect: float = 2.0,
    depth: int = 10_000,
    seed: int = 0,
) -> dict:
    rng = np.random.default_rng(seed)
    base = rng.lognormal(mean=2.0, sigma=1.0, size=n_taxa)
    base /= base.sum()

    group = rng.integers(0, 2, n_samples)
    counts = np.vstack([rng.multinomial(depth, base) for _ in range(n_samples)]).astype(float)

    planted_idx = rng.choice(n_taxa, planted, replace=False)
    counts[group == 1][:, planted_idx] *= effect

    taxa = [f"g__Taxon{i:03d}" for i in range(n_taxa)]
    samples = [f"S{i:03d}" for i in range(n_samples)]

    table = pd.DataFrame(counts, index=samples, columns=taxa, dtype=float)
    metadata = pd.DataFrame({"group": group.astype(str)}, index=samples)
    truth = {
        "planted_taxa": [taxa[i] for i in planted_idx],
        "effect": effect,
        "group_col": "group",
        "group_a": "0",
        "group_b": "1",
    }
    return {"counts": table, "metadata": metadata, "truth": truth}
