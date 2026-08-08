# Source-data inventory

This directory contains the numerical inputs used to construct the manuscript
figures. Files are grouped by their final figure assignment rather than by the
date on which an analysis was run.

## Main figures

- `main_figures/figure1/`: PC1, ATAC signal, gene density, GC content and loop
  tracks, plus the HiCExplorer track configuration. Whole-genome and regional
  contact matrices are excluded because of size and are regenerated from the
  deposited Hi-C reads.
- `main_figures/figure2/`: two-method and three-method trajectory calls,
  CLOUD assignments, class counts, percentages and pairwise tests.
- `main_figures/figure3/`: tissue-specificity and Ka/Ks values with plotting
  and statistical-test tables.
- `main_figures/figure4/`: TSS-centered ATAC profiles, per-gene chromatin
  features, compartment/TAD/loop summaries and tests.
- `main_figures/figures5_7/`: phylogenetic thresholds, Ks, gene age,
  expression, ATAC, compartment, TAD and contact results.
- `main_figures/figure8/`: regulatory-capture states, direction tests,
  adjusted models, TAD co-residence, coexpression and haplotype validation.

## Supplementary figures

`supplementary_figures/` contains the source and summary tables for
Supplementary Figures S1–S18. The accompanying `Supplementary_Sources/hapA/`
directory retains the replicate contact-decay, TAD-domain, TAD-boundary ATAC
and loop-associated gene inputs used for the replicate-level panels.

## Integrity manifest

`SOURCE_DATA_INDEX.tsv` lists every file under `source_data/`, its size and its
SHA-256 checksum. It can be regenerated after intentional source-data changes.
