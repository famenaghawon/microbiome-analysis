# Data Provenance

This document records exactly where the real dataset in this repository
came from, how it was processed, and how to reproduce it. Anyone can
follow these steps and obtain the same `data/cache/_genus_raw.tsv`.

## Source

- **BioProject:** NCBI SRA, real 16S rRNA-seq gut metagenome samples
- **Number of runs:** 124 (SRR41063392 – SRR41063515)
- **Organism:** human gut metagenome (TaxID 749906)
- **Platform:** Illumina MiSeq, paired-end
- **Sequencing strategy:** 16S rRNA gene amplicon

## Processing pipeline

The raw data was processed with **QIIME 2 (amplicon 2026.7)**:

1. Import paired-end reads (SRA FASTQ → QIIME 2 artifact)
2. Denoise with **DADA2** (standard 16S parameters)
3. Assign taxonomy with a **Greengenes2** reference
4. Collapse ASVs to genus rank → `all-genus-table.qza`
5. Export to BIOM → convert to TSV

The final commands (from the working directory of the QIIME 2 environment):

```bash
qiime tools export \
    --input-path all-genus-table.qza \
    --output-path genus-raw-export

biom convert \
    -i genus-raw-export/feature-table.biom \
    -o genus_raw.tsv \
    --to-tsv

cp genus_raw.tsv data/cache/_genus_raw.tsv
