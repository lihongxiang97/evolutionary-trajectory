# T2T pear 3D-genome analysis of duplicate genes

This repository contains the analysis, statistical outputs and figure source
data supporting **“Duplication-mode-specific 3D chromatin trajectories of
duplicate genes in pear.”** The gap-free, haplotype-resolved `Dangshansuli`
pear hapA assembly is the primary coordinate system; hapB is used to test
cross-haplotype robustness.

The study integrates Hi-C, ATAC-seq, fruit RNA-seq, multi-tissue expression,
comparative genomics and duplicate-gene evolutionary classifications.

## Analysis design

1. Process biological-replicate Hi-C, ATAC-seq and RNA-seq datasets against
   the gap-free hapA and hapB assemblies.
2. Identify A/B compartments, TADs, chromatin loops and OCRs using
   replicate-aware quality control.
3. Identify tandem, proximal and transposed duplicate pairs and polarize
   parent and child copies using peach synteny.
4. Define the primary four-class trajectory set by exact agreement between
   the expression-localization method and CDROM: `Con`, `Neo(C)`, `Neo(P)` and
   `Spe`.
5. All three trajectory algorithms can assign a subfunctionalization (`Sub`)
   state. `Sub` is absent from the primary two-method four-class set because
   the localization–CDROM intersection retains only exact matches among the
   four harmonized classes. The three-method sensitivity set likewise
   requires an exact four-class match across localization, CDROM and CLOUD;
   `Sub` and other discordant calls are classified as non-agreement rather
   than merged into another trajectory.
6. Test expression, accessibility, compartment, TAD, loop and
   promoter-distal OCR-capture associations using matched controls,
   covariate-adjusted models, effect sizes, confidence intervals and
   multiple-testing correction.
7. Quantify continuous-feature concordance, binary-state agreement and
   Young–Old effect-direction concordance through stringent reciprocal-best-
   hit hapA–hapB mappings. These comparisons assess robustness across two
   coordinate systems and are not treated as independent biological
   replication.

## Repository layout

```text
.
├── scripts/                              Analysis and figure-generation code
├── source_data/
│   ├── main_figures/                     Figure 1–8 plotting data and tracks
│   └── supplementary_figures/            Supplementary Figure S1–S18 sources
├── analysis_results/
│   ├── 3D_duplicate_mechanisms/          Pair- and gene-level 3D analyses
│   ├── Cross_Haplotype_Robustness/        hapA–hapB concordance analyses
│   └── Risk_Audit/                        TE/annotation sensitivity analyses
├── FIGURE_SOURCE_MAP.tsv                  Figure-to-data-to-script mapping
├── requirements.txt                       Python dependencies
└── README.md
```

`source_data/SOURCE_DATA_INDEX.tsv` records the path, file size and SHA-256
checksum for every figure-source file. Raw reads, genome assemblies, BAMs and
large Hi-C matrices are not stored in GitHub; their public accessions and
provenance are provided in the manuscript and Supplementary Table 1.

## Figure reproduction

- Figure 1: `scripts/make_figure1.py`
- Figure 2: `scripts/make_figure2_two_method.py`
- Figure 3: `scripts/prepare_figure3_data.py` and
  `scripts/make_figure3_two_method.py`
- Figure 4: `scripts/prepare_figure4_data.py` and
  `scripts/make_figure4_two_method.py`
- Figures 5–7: `scripts/make_figure5_revised_v2.py`,
  `scripts/make_figure6.py` and `scripts/make_figure7.py`
- Figure 8: `scripts/make_new_framework_figures.py`
- Supplementary Figure S1: `scripts/make_supplementary_figure_s1.py`
- Supplementary Figures S2–S18: the supplementary/risk/framework scripts
  listed in `FIGURE_SOURCE_MAP.tsv`

Each plotting script reads the corresponding directory under `source_data`
and writes regenerated files under `output/`. Figure 1 requires the large
merged and replicate Hi-C matrices listed in `FIGURE_SOURCE_MAP.tsv`; compact
tracks and loop tables are included here.

## Installation and checks

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts/build_source_index.py
python scripts/audit_repository_source_data.py
```

The upstream genomic workflow additionally uses HiC-Pro, HiCExplorer,
FitHiC2, standard alignment tools and interval-manipulation utilities.

## Interpretation safeguards

- Primary 3D calls use merged biological replicates; replicate-specific
  results are retained for concordance tests.
- Benjamini–Hochberg correction is applied within defined hypothesis
  families.
- Trans-loop results are reported descriptively because their replicate
  support is lower than that of cis loops.
- Promoter-distal OCR contacts are regulatory-capture signatures, not direct
  enhancer-validation experiments.
- hapB analyses test cross-coordinate-system robustness and do not constitute
  an independent biological experiment.

## Data availability

The 15-day-after-flowering pear Hi-C, ATAC-seq and RNA-seq datasets are
available from NGDC-GSA BioProject `PRJCA020047`, accession `CRA012754`.
Multi-tissue accessions, assembly references and dataset provenance are listed
in the manuscript and Supplementary Table 1.

## Citation

Please cite the accompanying manuscript and the published gap-free,
haplotype-resolved Asian and European pear genome assemblies when using this
workflow. The final manuscript citation will be added after publication.
