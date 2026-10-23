"""Headless-safe plotting helpers for microbiome analysis."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe: no display required
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def alpha_plot(alpha_df: pd.DataFrame, out_path: str) -> str:
    """Barplot of alpha diversity metrics per sample."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    metrics = ["richness", "shannon", "simpson"]
    titles = ["Richness", "Shannon", "Simpson"]
    colors = ["#2E86AB", "#A23B72", "#F18F01"]

    for ax, metric, title, color in zip(axes, metrics, titles, colors):
        ax.hist(alpha_df[metric], bins=20, color=color, edgecolor="white")
        ax.set_title(title)
        ax.set_xlabel(metric)
        ax.set_ylabel("count (samples)")

    fig.suptitle(f"Alpha diversity per sample (n={len(alpha_df)})")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def pcoa_plot(
    coords: pd.DataFrame,
    explained: np.ndarray,
    labels: pd.Series,
    out_path: str,
    title: str = "PCoA",
) -> str:
    """Scatter plot of PCoA coordinates, colored by group."""
    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    aligned = labels.loc[coords.index]
    unique = list(dict.fromkeys(aligned.values))
    cmap = plt.get_cmap("tab10")

    for i, group in enumerate(unique):
        mask = aligned == group
        ax.scatter(
            coords.loc[mask, "PC1"],
            coords.loc[mask, "PC2"],
            label=str(group),
            color=cmap(i % 10),
            alpha=0.75,
            edgecolor="k",
            linewidth=0.3,
            s=40,
        )

    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    ax.set_title(title)
    ax.legend(title="Group", frameon=False, fontsize=9)
    ax.axhline(0, color="grey", linewidth=0.4, linestyle=":")
    ax.axvline(0, color="grey", linewidth=0.4, linestyle=":")

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path
