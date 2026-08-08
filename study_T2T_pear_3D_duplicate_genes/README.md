# T2T pear 3D-genome analysis of duplicate genes

## Scope

This directory records the analysis workflow supporting the manuscript
*Duplication-mode-specific 3D chromatin trajectories of duplicate genes in
pear*. The gap-free `Dangshansuli` hapA assembly is the discovery coordinate
system and hapB is used as an independent haplotype-resolved validation
coordinate system.

The workflow integrates Hi-C, ATAC-seq, fruit RNA-seq, multi-tissue expression,
comparative genomics and duplicate-gene classifications. The deposited files
are analysis scripts and compact source tables; raw sequencing data remain in
the public accessions reported in the manuscript.

## Analysis logic

1. Align and process Hi-C, ATAC-seq and RNA-seq against the gap-free hapA and
   hapB assemblies.
2. Call compartments, TADs, loops and open chromatin using replicate-aware
   workflows.
3. Identify tandem, proximal and transposed duplicate pairs and polarize parent
   and child copies using synteny with peach.
4. Define the primary trajectory set by exact agreement between the
   expression-localization method and CDROM on four labels: `Con`, `Neo(C)`,
   `Neo(P)` and `Spe`.
5. Use CLOUD as an independent sensitivity classifier. CLOUD's `Sub` label is
   retained for posterior-confidence assessment but is not merged with any of
   the four harmonized classes. A pair enters the strict three-method set only
   when CLOUD assigns the same four-class label; `Sub` and other discordant
   calls are counted as non-agreement.
6. Quantify expression, accessibility, compartment, TAD, loop and
   promoter-distal OCR-capture associations using matched controls,
   covariate-adjusted models and multiple-testing correction.
7. Validate continuous features, binary 3D states and Young-Old effect
   directions in stringent reciprocal-best-hit hapA-hapB mappings.

## Contact-decay display rule

Supplementary Figure S1 excludes the low-support extreme-distance tail. For
each map independently, only distance bins supported by at least 1,000 possible
chromosome-bin pairs and at least 10% of that map's maximum possible-pair count
are plotted. The threshold is applied to support counts, not selected from the
shape of the resulting curve.

## Directory contents

- `scripts/`: analysis, figure-generation and package-audit scripts.
- `source_data/figure2/`: four-class consensus counts, percentages, statistical
  tests and plot data.
- `source_data/supplementary_figure_s1/`: filtered contact-decay data and
  replicate-to-merged correlations.
- `source_data/hapB_validation/`: the seven result blocks reported in
  Supplementary Table 17.

## Key software

Hi-C processing used HiC-Pro, HiCExplorer and FitHiC2. Downstream analysis and
figures use Python with NumPy, pandas, SciPy, statsmodels, matplotlib and
seaborn. Exact versions should be recorded from the analysis environments used
for the final run.

## Reproducibility notes

- Merged biological replicates are used for primary feature calling and
  visualization; replicate-specific results are retained for concordance and
  robustness analyses.
- Benjamini-Hochberg correction is applied where families of hypotheses are
  tested.
- Trans-loop findings are treated as descriptive because replicate support is
  substantially lower than for cis loops.
- Promoter-distal OCR contacts are regulatory-capture signatures, not direct
  enhancer-validation assays.

## Data access

The 15-day-after-flowering pear Hi-C, ATAC-seq and RNA-seq data are available
from NGDC-GSA BioProject `PRJCA020047`, accession `CRA012754`. Multi-tissue
accessions and assembly references are listed in the manuscript and
Supplementary Table 1.
