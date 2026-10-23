"""Microbiome Analysis Dashboard.

Interactive Streamlit app that loads real HMP-adjacent data (or user uploads)
and runs live analysis: alpha/beta diversity, PCoA, DA, and ML classification.

Design principle: this is a *view layer* over the same modules used in
scripts/run_pipeline.py and scripts/run_ml.py. No logic duplicates.
"""
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from microbiome_analysis.data.loaders import load
from microbiome_analysis.differential.wilcoxon import wilcoxon_da
from microbiome_analysis.diversity.alpha import alpha_table
from microbiome_analysis.diversity.beta import beta_distances, distances_to_frame
from microbiome_analysis.diversity.ordination import ordination_coords
from microbiome_analysis.diversity.permanova import run_permanova
from microbiome_analysis.ml.classifier import evaluate, train_final_and_rank_features
from microbiome_analysis.preprocess.transforms import clr

# ─────────────────────────────────────────────────────────────────────────────
# Page config + state
# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Microbiome Analysis",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .block-container { padding-top: 2rem; }
    .metric-big { font-size: 2.2rem; font-weight: 700; }
    .metric-label { font-size: 0.85rem; opacity: 0.7; text-transform: uppercase; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    </style>
    """,
    unsafe_allow_html=True,
)


# ─────────────────────────────────────────────────────────────────────────────
# Data cache (avoid reloading on every widget interaction)
# ─────────────────────────────────────────────────────────────────────────────
@st.cache_data(show_spinner="Loading microbiome data...")
def load_data(source: str, max_samples: int | None):
    bundle = load(source, path=None, max_samples=max_samples)
    return (
        bundle["counts"],
        bundle["metadata"],
        bundle.get("report", {}),
        bundle.get("source", source),
    )


@st.cache_data(show_spinner="Computing alpha diversity...")
def compute_alpha(counts_hash: str, counts: pd.DataFrame):
    return alpha_table(counts)


@st.cache_data(show_spinner="Computing distances...")
def compute_beta(counts_hash: str, counts: pd.DataFrame):
    return beta_distances(counts)


@st.cache_data(show_spinner="Running ordination...")
def compute_ordination(counts_hash: str, counts: pd.DataFrame, metric: str):
    dists = beta_distances(counts)
    dm = dists[metric]
    coords, explained = ordination_coords(dm)
    return coords, explained


@st.cache_data(show_spinner="Running PERMANOVA...")
def compute_permanova(counts_hash: str, counts: pd.DataFrame,
                      labels_tuple, metric: str):
    dists = beta_distances(counts)
    dm = dists[metric]
    labels = pd.Series(dict(labels_tuple))
    return run_permanova(dm, labels, permutations=999)


@st.cache_data(show_spinner="Running differential abundance...")
def compute_da(counts_hash: str, counts: pd.DataFrame,
               labels_tuple, ga: str, gb: str, method: str):
    labels = pd.Series(dict(labels_tuple))
    if method == "wilcoxon":
        return wilcoxon_da(clr(counts), labels, group_a=ga, group_b=gb)
    else:
        from microbiome_analysis.differential.ancom_bc import ancom_bc_da
        return ancom_bc_da(counts, labels, ref_group=ga)


@st.cache_data(show_spinner="Training leakage-safe classifier...")
def compute_ml(counts_hash: str, counts: pd.DataFrame,
               labels_tuple, k: int, folds: int):
    labels = pd.Series(dict(labels_tuple))
    metrics = evaluate(counts, labels, k=k, folds=folds, seed=42)
    importance = train_final_and_rank_features(counts, labels, k=min(k, 30), seed=42)
    return metrics, importance


# ─────────────────────────────────────────────────────────────────────────────
# Sidebar - data source
# ─────────────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🦠 Microbiome Analysis")
    st.caption("Reproducible analysis of real HMP-adjacent data")

    st.subheader("Data source")
    source = st.selectbox(
        "Load dataset",
        [
            "gmrepo_v2 (new, multi-site, 662 samples)",
            "gmrepo (old, 122 samples)",
            "synthetic (test fixture)",
        ],
        index=0,
    )
    if source.startswith("gmrepo_v2"):
        source_key = "gmrepo_v2"
    elif source.startswith("gmrepo"):
        source_key = "gmrepo"
    else:
        source_key = "synthetic"
    max_samples = st.slider(
        "Max samples (for responsiveness)",
        min_value=20, max_value=200, value=122, step=2,
    )

    st.divider()

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────
counts, meta, report, source_label = load_data(source_key, max_samples)

# Sanity: align counts and meta
meta = meta.reindex(counts.index)
meta = meta.fillna("missing")

# Compute depth OR relative abundance metric
depths = counts.sum(axis=1)
med = depths.median()

# Detect if data is already relative (row sums ~ 1.0)
is_relative = med < 10
meta["depth_bin"] = ["high" if d >= med else "low" for d in depths]
PREFERRED_GROUP_COLS = ["BodySiteCategory", "BodySite", "group", "depth_bin"]

# Counts hash for caching
counts_hash = f"{counts.shape[0]}_{counts.shape[1]}_{int(depths.sum())}"

# ─────────────────────────────────────────────────────────────────────────────
# Sidebar - grouping controls
# ─────────────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.subheader("Analysis controls")

    # Only offer columns with 2-10 unique values (categorical-friendly)
    categorical_cols = []
    for c in meta.columns:
        n_unique = meta[c].nunique()
        if 2 <= n_unique <= 10:
            categorical_cols.append(c)

    # Order by preference (real biology first, depth_bin last)
    preferred_order = ["BodySiteCategory", "BodySite", "group", "depth_bin"]
    categorical_cols.sort(
        key=lambda c: preferred_order.index(c) if c in preferred_order else 999
    )

    if not categorical_cols:
        st.error("No metadata column has ≥2 groups. Defaulting to depth_bin.")
        group_col = "depth_bin"
    else:
        group_col = st.selectbox(
            "Group samples by",
            categorical_cols,
            index=0,
        )

    # Show warning if chosen column has only 1 unique value
    n_groups = meta[group_col].nunique()
    if n_groups < 2:
        st.warning(
            f"'{group_col}' has only {n_groups} unique value(s). "
            "PERMANOVA and DA need ≥ 2 groups."
        )

    # Group A / B selectors only if 2+ groups
    if n_groups >= 2:
        unique_groups = sorted(meta[group_col].astype(str).unique())
        
        # Auto-detect the best default contrast
        default_a_idx = 0
        default_b_idx = 1 if len(unique_groups) > 1 else 0
        
        # If BodySiteCategory is selected, automatically default to Gut vs Oral
        if "Gut" in unique_groups and "Oral" in unique_groups:
            default_a_idx = unique_groups.index("Gut")
            default_b_idx = unique_groups.index("Oral")
            
        group_a = st.selectbox("Group A", unique_groups, index=default_a_idx)
        group_b = st.selectbox("Group B", unique_groups, index=default_b_idx)

        # Helpful hint when using body site
        if group_col == "BodySiteCategory" and "Gut" in unique_groups and "Oral" in unique_groups:
            st.caption(
                "💡 **Tip:** `Gut` vs `Oral` is the classic microbiome contrast — "
                "well-powered for ML and DA. The groups are pre-selected for you."
            )
    else:
        group_a = group_b = None

    st.divider()
# ─────────────────────────────────────────────────────────────────────────────
# Header
# ─────────────────────────────────────────────────────────────────────────────
st.title("Microbiome Analysis — Reproducible by Design")
st.caption(
    f"Source: **{source_label}** · {counts.shape[0]} samples × "
    f"{counts.shape[1]} genera · median depth {depths.median():,.0f} reads"
)

# Big metrics row
c1, c2, c3, c4 = st.columns(4)
with c1:
    st.markdown(
        f'<div class="metric-big">{counts.shape[0]:,}</div>'
        f'<div class="metric-label">Samples</div>',
        unsafe_allow_html=True,
    )
with c2:
    st.markdown(
        f'<div class="metric-big">{counts.shape[1]:,}</div>'
        f'<div class="metric-label">Genera</div>',
        unsafe_allow_html=True,
    )
with c3:
    # Show mean relative abundance if data is normalized, else depth
    if is_relative:
        val = f"{depths.mean():.2f}"
        label = "Mean Rel. Abund."
    else:
        val = f"{depths.median():,.0f}"
        label = "Median depth"
    st.markdown(
        f'<div class="metric-big">{val}</div>'
        f'<div class="metric-label">{label}</div>',
        unsafe_allow_html=True,
    )
with c4:
    if is_relative:
        val = f"{counts.shape[1]}"
        label = "Total Genera"
    else:
        val = f"{depths.sum():,.0f}"
        label = "Total reads"
    st.markdown(
        f'<div class="metric-big">{val}</div>'
        f'<div class="metric-label">{label}</div>',
        unsafe_allow_html=True,
    )

st.divider()

# ─────────────────────────────────────────────────────────────────────────────
# Tabs
# ─────────────────────────────────────────────────────────────────────────────
tab_overview, tab_alpha, tab_beta, tab_da, tab_ml, tab_data = st.tabs([
    "📊 Overview",
    "📈 Alpha diversity",
    "🌐 Beta diversity",
    "🔬 Differential abundance",
    "🤖 Machine learning",
    "📋 Data & QC",
])

# ─── Overview ────────────────────────────────────────────────────────────────
with tab_overview:
    st.subheader("Top genera by total reads")
    totals = counts.sum().sort_values(ascending=False)
    top20 = totals.head(20)
    fig = px.bar(
        x=top20.values,
        y=top20.index,
        orientation="h",
        labels={"x": "Total reads", "y": "Genus"},
        title=f"Top 20 genera (out of {counts.shape[1]})",
        color=top20.values,
        color_continuous_scale="Viridis",
    )
    fig.update_layout(
        yaxis={"categoryorder": "total ascending"},
        showlegend=False,
        coloraxis_showscale=False,
        height=550,
    )
    st.plotly_chart(fig, width="stretch")

    st.subheader("Sample depth distribution")
    fig = px.histogram(
        x=depths.values,
        nbins=30,
        labels={"x": "Sequencing depth (reads)"},
        title=f"Depth distribution (n={len(depths)}, median={depths.median():,.0f})",
    )
    fig.update_traces(marker_color="#2E86AB")
    st.plotly_chart(fig, width="stretch")

# ─── Alpha diversity ─────────────────────────────────────────────────────────
with tab_alpha:
    st.subheader("Alpha diversity per sample")
    alpha = compute_alpha(counts_hash, counts)
    alpha["sample"] = alpha.index
    alpha["group"] = meta[group_col].values

    metric = st.radio(
        "Metric",
        ["richness", "shannon", "simpson"],
        horizontal=True,
    )

    fig = px.box(
        alpha,
        x="group",
        y=metric,
        color="group",
        points="all",
        title=f"{metric.title()} by {group_col}",
        labels={metric: metric.title(), "group": group_col},
    )
    fig.update_layout(showlegend=False, height=450)
    st.plotly_chart(fig, width="stretch")

    c1, c2 = st.columns(2)
    with c1:
        st.metric("Mean", f"{alpha[metric].mean():.3f}")
    with c2:
        st.metric("Std", f"{alpha[metric].std():.3f}")

    with st.expander("View raw alpha table"):
        st.dataframe(alpha.set_index("sample"), width="stretch")

# ─── Beta diversity ──────────────────────────────────────────────────────────
with tab_beta:
    st.subheader("Beta diversity ordination")

    if n_groups < 2:
        st.error(
            f"PERMANOVA requires ≥ 2 groups in '{group_col}'. "
            "Pick a different grouping column in the sidebar."
        )
    else:
        metric = st.radio("Distance metric", ["braycurtis", "jaccard"], horizontal=False)

        coords, explained = compute_ordination(counts_hash, counts, metric)
        coords["sample"] = coords.index
        coords["group"] = meta[group_col].values

        fig = px.scatter(
            coords,
            x="PC1",
            y="PC2",
            color="group",
            hover_data={"sample": True, "PC1": ":.3f", "PC2": ":.3f"},
            title=f"PCoA — {metric.title()} "
                  f"(PC1={explained[0]*100:.1f}%, PC2={explained[1]*100:.1f}%)",
            color_discrete_sequence=px.colors.qualitative.Set2,
        )
        fig.update_traces(marker={"size": 10, "line": {"width": 1, "color": "white"}})
        fig.update_layout(height=550)
        st.plotly_chart(fig, width="stretch")

        # PERMANOVA
        labels_tuple = tuple(zip(meta.index, meta[group_col].astype(str)))
        permanova_result = compute_permanova(counts_hash, counts, labels_tuple, metric)

        c1, c2, c3 = st.columns(3)
        c1.metric("Pseudo-F", f"{permanova_result['pseudo_f']:.2f}")
        c2.metric("p-value", f"{permanova_result['p_value']:.4f}")
        c3.metric("Permutations", f"{permanova_result['permutations']}")

        if permanova_result["p_value"] < 0.05:
            st.success(
                f"PERMANOVA: communities differ significantly by `{group_col}` "
                f"(p = {permanova_result['p_value']:.4f})"
            )
        else:
            st.info(
                f"PERMANOVA: no significant difference by `{group_col}` "
                f"(p = {permanova_result['p_value']:.4f})"
            )

        with st.expander("Download distance matrix"):
            dm = distances_to_frame(beta_distances(counts)[metric])
            st.download_button(
                "Download DM (TSV)",
                dm.to_csv(sep="\t").encode(),
                file_name=f"distance_matrix_{metric}.tsv",
                mime="text/tab-separated-values",
            )

# ─── Differential abundance ──────────────────────────────────────────────────
with tab_da:
    st.subheader("Differential abundance")

    # Count samples per group to ensure DA is valid
    group_counts = meta[group_col].astype(str).value_counts()
    n_a = int(group_counts.get(str(group_a), 0))
    n_b = int(group_counts.get(str(group_b), 0))

    if group_a == group_b:
        st.error("Group A and Group B must be different.")
    elif n_a < 3 or n_b < 3:
        st.error(
            f"DA requires ≥3 samples per group. "
            f"Found: {group_a} (n={n_a}), {group_b} (n={n_b}). "
            f"Please choose larger groups (e.g., Gut vs Oral)."
        )
    else:
        c1, c2 = st.columns(2)
        with c1:
            method = st.selectbox("Method", ["wilcoxon", "ancom_bc"])
        with c2:
            alpha_thresh = st.slider("FDR threshold (q <)", 0.01, 0.20, 0.05, 0.01)

        labels_tuple = tuple(zip(meta.index, meta[group_col].astype(str)))
        da = compute_da(counts_hash, counts, labels_tuple, group_a, group_b, method)

        n_sig = (da["q"] < alpha_thresh).sum()
        c1, c2, c3 = st.columns(3)
        c1.metric("Significant", f"{n_sig} / {len(da)}")
        c2.metric("Group A", f"{group_a} (n={n_a})")
        c3.metric("Group B", f"{group_b} (n={n_b})")

        # Volcano plot
        da_plot = da.copy()
        da_plot["neg_log10_q"] = -np.log10(da_plot["q"].clip(lower=1e-300))
        da_plot["significant"] = da_plot["q"] < alpha_thresh
        da_plot["direction"] = np.where(da_plot["effect"] > 0, "up in B", "up in A")
        da_plot.loc[~da_plot["significant"], "direction"] = "not significant"

        fig = px.scatter(
            da_plot,
            x="effect",
            y="neg_log10_q",
            color="direction",
            hover_data={"taxon": True, "effect": ":.3f", "q": ":.2e"},
            title=f"Volcano plot — {method} ({group_a} vs {group_b})",
            color_discrete_map={
                "up in A": "#e74c3c",
                "up in B": "#3498db",
                "not significant": "#cccccc",
            },
        )
        fig.add_hline(y=-np.log10(alpha_thresh), line_dash="dash", line_color="grey")
        fig.add_vline(x=0, line_dash="dot", line_color="grey")
        fig.update_layout(height=500)
        st.plotly_chart(fig, width="stretch")

        st.subheader(f"Top {min(20, n_sig if n_sig else len(da))} by q-value")
        st.dataframe(
            da.head(20)[["taxon", "effect", "p", "q", "significant"]],
            width="stretch",
        )

        st.download_button(
            "Download full DA table (TSV)",
            da.to_csv(sep="\t", index=False).encode(),
            file_name=f"da_{method}_{group_a}_vs_{group_b}.tsv",
            mime="text/tab-separated-values",
        )

# ─── ML ──────────────────────────────────────────────────────────────────────
with tab_ml:
    st.subheader("Leakage-safe classification")

    # Count samples per group
    group_counts = meta[group_col].astype(str).value_counts()
    min_group_size = int(group_counts.min())
    n_groups_for_ml = len(group_counts)

    if n_groups_for_ml < 2:
        st.error(
            f"ML requires ≥ 2 groups in '{group_col}'. "
            "Pick a different grouping column in the sidebar."
        )
    else:
        st.caption(
            f"Classes: {dict(group_counts)} · "
            f"smallest = {min_group_size} samples"
        )

        # If any group has < 2 samples, we cannot run ML at all
        if min_group_size < 2:
            st.error(
                f"ML cannot run: '{group_counts.idxmin()}' has only "
                f"{min_group_size} sample(s). "
                f"Switch to `Group A: Gut` / `Group B: Oral` in the sidebar."
            )
        else:
            c1, c2 = st.columns(2)
            with c1:
                k_features = st.slider(
                    "Top-k features", 10, min(100, counts.shape[1]), 50, 5
                )
            with c2:
                # Safely compute slider bounds
                max_folds = max(2, min(10, min_group_size))
                min_folds = 2
                
                if min_folds >= max_folds:
                    st.info(f"Only {min_group_size} samples in smallest class. Using 2-fold CV.")
                    folds = 2
                else:
                    default_folds = min(5, max_folds)
                    folds = st.slider("CV folds", min_folds, max_folds, default_folds)

            labels_tuple = tuple(zip(meta.index, meta[group_col].astype(str)))
            try:
                metrics, importance = compute_ml(
                    counts_hash, counts, labels_tuple, k_features, folds
                )

                c1, c2, c3, c4 = st.columns(4)
                c1.metric(
                    "Balanced accuracy",
                    f"{metrics['balanced_accuracy_mean']:.3f}",
                    delta=f"±{metrics['balanced_accuracy_std']:.3f}",
                )
                if "auc_mean" in metrics:
                    c2.metric(
                        "AUC",
                        f"{metrics['auc_mean']:.3f}",
                        delta=f"±{metrics['auc_std']:.3f}",
                    )
                c3.metric("Classes", metrics["n_classes"])
                c4.metric("Folds", folds)

                ba = metrics["balanced_accuracy_mean"]
                if ba > 0.95:
                    st.error("⚠️ Balanced accuracy > 0.95 — check for leakage.")
                elif ba > 0.6:
                    st.success(
                        f"✓ Balanced accuracy {ba:.3f} is in the expected range."
                    )
                else:
                    st.info(f"ℹ️ Balanced accuracy {ba:.3f} is near chance (0.5).")

                st.subheader("Feature importance (top 20)")
                top_imp = importance.head(20).sort_values("importance")
                fig = px.bar(
                    top_imp,
                    x="importance",
                    y="taxon",
                    orientation="h",
                    title="Which genera matter most for prediction?",
                    color="importance",
                    color_continuous_scale="Teal",
                )
                fig.update_layout(
                    showlegend=False, coloraxis_showscale=False, height=550
                )
                st.plotly_chart(fig, width="stretch")

                with st.expander("View full feature importance table"):
                    st.dataframe(importance, width="stretch")
            except Exception as e:  # noqa: BLE001
                st.error(
                    f"ML failed with {folds}-fold CV: {type(e).__name__}: {e}\n\n"
                    f"Try: lower folds to {max(2, min_group_size - 1)}, or pick a 2-class "
                    f"contrast (Gut vs Oral recommended)."
                )

# ─── Data & QC ───────────────────────────────────────────────────────────────
with tab_data:
    st.subheader("Loader QC report")
    if report:
        st.json(report)
    else:
        st.info("No QC report.")

    st.subheader("Metadata preview")
    st.dataframe(meta.head(20), width="stretch")

    st.subheader("Counts preview (5 × 10)")
    st.dataframe(counts.iloc[:5, :10], width="stretch")

    st.subheader("Download raw counts")
    st.download_button(
        "Download counts (TSV)",
        counts.to_csv(sep="\t").encode(),
        file_name="counts.tsv",
        mime="text/tab-separated-values",
    )

st.caption(
    "Built with Streamlit · Plotly · scikit-bio · "
    "Source: github.com/famenaghawon/microbiome-analysis"
)
