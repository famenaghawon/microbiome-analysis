# Data format contract

## `hmp_genus.tsv`
- Tab-separated, first column = sample ID (used as index).
- Remaining columns = genus names prefixed `g__` (e.g. `g__Bacteroides`).
- Values = integer read counts.

## `hmp_metadata.tsv`
- Tab-separated, first column = sample ID (must intersect table index ≥ 80%).
- Required column `body_site` ∈ {gut, oral, skin, nasal, vaginal}.
  Loader lowercases and strips.

## Suggested sources (all public, no auth)
- **HMP1-II 16S Qiita study 1928** — export genus-level table + metadata via Qiita UI.
- **`curatedMetagenomicData` (Bioconductor)** — `HMP_2012.genus` (WGS).
- **`microbiome` R package, `hmp2`** — genus-level, small.

Drop the TSVs into `data/cache/` (gitignored) or host them and set `$env:MB_HMP_URL`.
