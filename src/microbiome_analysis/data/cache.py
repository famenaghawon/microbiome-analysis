# data/cache.py
CACHE = Path(os.getenv("MB_CACHE", "data/cache"))

def cached(name, url=None, loader=None, force=False):
    p = CACHE / name
    if p.exists() and not force:
        return p
    p.parent.mkdir(parents=True, exist_ok=True)
    if url:
        urllib.request.urlretrieve(url, p)
    elif loader:
        loader().to_csv(p, sep="\t", index=False)
    return p

# data/loaders.py
def load_hmp_subset():
    """~150 samples, gut/oral/skin genera — zenodo-hosted, no auth."""
    tsv = cached("hmp_subset_counts.tsv",
                 url="https://zenodo.org/records/XXXX/files/hmp_subset.tsv")
    return pd.read_csv(tsv, sep="\t", index_col=0)

def validate(table, metadata):
    inter = table.index.intersection(metadata.index)
    report = {
        "n_table": len(table), "n_meta": len(metadata), "n_joined": len(inter),
        "join_rate": len(inter) / len(table),
        "median_depth": float(table.loc[inter].sum(axis=1).median()),
        "taxa_with_zero_prevalence": int((table.loc[inter].sum(0) == 0).sum()),
    }
    if report["join_rate"] < 0.8:
        raise ValueError(f"Join rate {report['join_rate']:.2f} < 0.8 — check IDs")
    return report

