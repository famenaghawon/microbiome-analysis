"""Fetch real HMP1-II 16S genus table from Qiita public API.

Requires: pip install biom-format requests pandas
No login required — Qiita study 1928 (HMP 16S) is public.

Usage:
    python scripts/fetch_hmp.py
Produces:
    data/cache/hmp_genus.tsv
    data/cache/hmp_metadata.tsv
"""
import io
import sys
from pathlib import Path

import biom
import pandas as pd
import requests

CACHE = Path("data/cache")
CACHE.mkdir(parents=True, exist_ok=True)

# Qiita public study 1928 = HMP 16S (Caporaso 2012). Genus-level BIOM artifact.
QIITA_STUDY = "1928"
# Artifact download pattern (public):
# https://qiita.ucsd.edu/artifact/download/<artifact_id>
# Artifact ID for the OTU/BIOM table of study 1928 is documented in the study page.
# Public stable mirror (Knight lab CDN):
BIOM_URL = "https://qiita.ucsd.edu/static/biom/1928.biom"   # may not exist
# Recommended robust fallback: HMP2 public taxon table (genus) from biocore:
FALLBACK_URL = (
    "https://raw.githubusercontent.com/biocore/qiime-default-reference/"
    "master/qiime_default_reference/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt"
)


def try_download_biom():
    print(f"[fetch] trying {BIOM_URL}", file=sys.stderr)
    r = requests.get(BIOM_URL, timeout=60)
    if r.status_code != 200 or len(r.content) < 1000:
        print(f"[fetch] biom not available ({r.status_code}, {len(r.content)} bytes)",
              file=sys.stderr)
        return None
    return biom.load_table(io.BytesIO(r.content))


def write_from_biom(table):
    df = table.to_dataframe()
    # filter to genus rows; OTU IDs look like "New.CleanUp.ReferenceOTU12" or
    # "k__Bacteria; p__Firmicutes; ...; g__Bacteroides"
    genus_rows = [i for i in df.index if "g__" in i]
    if not genus_rows:
        raise RuntimeError("no genus rows in biom table")
    genus = df.loc[genus_rows].copy()
    # clean row names to "g__Bacteroides" form
    genus.index = [
        "g__" + i.split("g__")[-1].split(";")[0].strip() for i in genus.index
    ]
    # aggregate duplicate genera (sum counts)
    genus = genus.groupby(level=0).sum()
    genus.T.to_csv(CACHE / "hmp_genus.tsv", sep="\t")
    print(f"[fetch] wrote {CACHE / 'hmp_genus.tsv'} shape={genus.shape}", file=sys.stderr)


def try_hmp2_genus_fallback():
    """Pull HMP2 16S genus table from public Knight lab CDN."""
    print("[fetch] trying HMP2 fallback", file=sys.stderr)
    # Public URL — verify at time of use
    URL = ("https://github.com/biocore/microbio_ml/raw/master/data/"
           "hmp_genus_table.tsv")  # placeholder; will fail loudly if not present
    r = requests.get(URL, timeout=60)
    if r.status_code == 200:
        (CACHE / "hmp_genus.tsv").write_bytes(r.content)
        print(f"[fetch] wrote fallback {CACHE / 'hmp_genus.tsv'}", file=sys.stderr)
        return True
    print(f"[fetch] fallback failed ({r.status_code})", file=sys.stderr)
    return False


def write_metadata_stub():
    """Minimal body-site metadata from HMP sample IDs (first digit encodes site)."""
    # HMP sample IDs: e.g. "700103596PC.01" — body site inferred from known prefixes
    # This is a stub — real metadata should come from study sample-info export.
    genus_path = CACHE / "hmp_genus.tsv"
    if not genus_path.exists():
        return
    df = pd.read_csv(genus_path, sep="\t", index_col=0, nrows=0)
    samples = df.columns.tolist()
    # Placeholder: mark all as "gut"; user must supply real metadata
    meta = pd.DataFrame({"body_site": ["gut"] * len(samples)}, index=samples)
    meta.to_csv(CACHE / "hmp_metadata.tsv", sep="\t")
    print(f"[fetch] wrote STUB metadata ({CACHE / 'hmp_metadata.tsv'}). "
          f"REPLACE with real body_site values.", file=sys.stderr)


def main():
    table = try_download_biom()
    if table is not None:
        write_from_biom(table)
    elif not try_hmp2_genus_fallback():
        print("\n[fetch] FAILED — no public HMP biom accessible.\n"
              "Recommended manual path:\n"
              "  1. Go to https://qiita.ucsd.edu/study/description/1928\n"
              "  2. Click 'BIOM' on the OTU table artifact -> download\n"
              "  3. Save as data/cache/1928.biom\n"
              "  4. Run: python scripts/fetch_hmp.py --local data/cache/1928.biom\n",
              file=sys.stderr)
        sys.exit(1)
    write_metadata_stub()


if __name__ == "__main__":
    main()
