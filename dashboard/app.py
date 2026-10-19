# dashboard/app.py
import streamlit as st
from microbiome_analysis.data.loaders import load_hmp_subset
from microbiome_analysis.diversity.alpha import alpha_table
from microbiome_analysis.pipeline import run_beta, run_da, run_ml

st.title("Microbiome analysis — reproducible by design")
src = st.sidebar.selectbox("Dataset", ["synthetic", "hmp_subset", "local BIOM"])
bundle = load_bundle(src)  # {counts, meta, tree?}

tab_a, tab_b, tab_d, tab_m = st.tabs(
    ["Alpha", "Beta / PCoA", "Differential abundance", "ML"])
with tab_a: st.dataframe(alpha_table(bundle["counts"]))
with tab_b: st.pyplot(run_beta(bundle["counts"], bundle["meta"]))
with tab_d: st.dataframe(run_da(bundle["counts"], bundle["meta"]))
with tab_m: st.json(run_ml(bundle["counts"], bundle["meta"]))
