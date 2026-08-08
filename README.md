# T2T pear 3D-genome analysis of duplicate genes

This repository contains the analysis and figure-generation code supporting
the study **“Duplication-mode-specific 3D chromatin trajectories of duplicate
genes in pear.”** It supersedes the earlier exploratory repository content.

The gap-free, haplotype-resolved `Dangshansuli` pear hapA assembly is used as
the discovery coordinate system, and hapB provides an independent
cross-haplotype validation system. The project integrates Hi-C, ATAC-seq,
fruit RNA-seq, multi-tissue expression, comparative genomics and
duplicate-gene evolutionary classifications.

## Scientific workflow

1. Process replicate Hi-C, ATAC-seq and RNA-seq datasets against the gap-free
   hapA and hapB assemblies.
2. Identify A/B compartments, TADs, chromatin loops and open chromatin regions
   with replicate-aware quality control.
3. Identify tandem, proximal and transposed duplicate pairs and polarize
   parent and child copies using peach synteny.
4. Define the primary four-class trajectory set by exact agreement between
   the expression-localization method and CDROM: `Con`, `Neo(C)`, `Neo(P)` and
   `Spe`.
5. Use CLOUD as an independent sensitivity classifier. CLOUD's
   subfunctionalization (`Sub`) state is retained for posterior-confidence
   assessment but is not merged into the four harmonized classes. The strict
   three-method set requires an exact four-class match; `Sub` and all other
   discordant calls are treated as non-agreement.
6. Test expression, accessibility, compartment, TAD, loop and promoter-distal
   OCR-capture associations using matched controls, covariate-adjusted models,
   effect sizes, confidence intervals and multiple-testing correction.
7. Validate quantitative features, binary 3D states and Young–Old effect
   directions using stringent reciprocal-best-hit hapA–hapB mappings.

## Repository layout

```text
.
├── scripts/                         Analysis, plotting and package-audit code
├── source_data/
│   ├── figure2/                     Four-class consensus source tables
│   ├── hapB_validation/             Supplementary Table 17 source results
│   └── supplementary_figure_s1/     Contact-decay source data and summary
├── requirements.txt                 Python dependencies
└── README.md
```

The repository intentionally contains compact source tables and executable
analysis code, not raw sequencing reads, genome assemblies or large Hi-C
matrices.

## Main scripts

- `scripts/analyze_haplotype_validation.py`: hapA–hapB feature and age-effect
  validation.
- `scripts/make_figure2_two_method.py`: Figure 2 using the primary four-class
  two-method consensus.
- `scripts/make_figure4_two_method.py`: duplicate-copy chromatin-state plots.
- `scripts/make_figure5_revised_v2.py`, `make_figure6.py`, `make_figure7.py`:
  downstream evolutionary and chromatin analyses.
- `scripts/make_supplementary_figures_v3.py`: replicate and robustness
  supplementary figures, including support-filtered contact decay.
- `scripts/make_risk_figures.py`, `make_te_sensitivity_figure.py` and
  `make_new_framework_figures.py`: matched-control, annotation, TE and
  mechanistic sensitivity analyses.
- `scripts/audit_submission_package.py` and
  `audit_plant_communications_package.py`: submission-package consistency
  checks.

## Contact-decay support rule

Supplementary Figure S1 excludes the low-support extreme-distance tail. For
each map independently, a distance bin is displayed only if it contains at
least 1,000 possible chromosome-bin pairs and at least 10% of that map's
maximum possible-pair count. This is a support-based rule defined independently
of curve shape.

## Installation

Create an isolated Python environment and install the listed dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The upstream genomic workflow additionally uses HiC-Pro, HiCExplorer,
FitHiC2, standard alignment tools and interval-manipulation utilities. These
large external programs and their reference indices are not vendored here.

## Reproducibility and interpretation

- Primary 3D feature calls and visualizations use merged biological
  replicates; replicate-specific results are retained for concordance tests.
- Benjamini–Hochberg correction is applied within defined hypothesis families.
- Trans-loop results are reported descriptively because their replicate
  support is lower than that of cis loops.
- Promoter-distal OCR contacts are interpreted as regulatory-capture
  signatures rather than direct enhancer-validation experiments.
- File paths and sample metadata are supplied through script arguments or
  configuration variables and should be adapted to the local installation.

## Data availability

The 15-day-after-flowering pear Hi-C, ATAC-seq and RNA-seq datasets are
available from NGDC-GSA BioProject `PRJCA020047`, accession `CRA012754`.
Multi-tissue accessions, assembly references and complete dataset provenance
are provided in the manuscript and Supplementary Table 1.

## Citation

Please cite the accompanying manuscript and the published gap-free,
haplotype-resolved Asian and European pear genome assemblies when using this
workflow. The final bibliographic citation will be added after publication.
