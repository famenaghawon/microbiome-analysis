"""Microbiome Analysis Dashboard.

Placeholder Streamlit app. Phase 4 will add real interactive tabs for
alpha diversity, beta ordination, differential abundance, and ML metrics.
"""
import streamlit as st

st.set_page_config(page_title="Microbiome Analysis", layout="wide")
st.title("Microbiome Analysis — Reproducible by Design")
st.markdown(
    """
    **Status:** Phase 2 complete (real-data DA pipeline + CI).
    Dashboard UI is Phase 4.

    For the current analysis, run:

    ```powershell
    python scripts/run_pipeline.py --source gmrepo --out results/real_sanity/
    ```
    """
)
