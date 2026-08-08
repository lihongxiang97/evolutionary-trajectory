# Evolutionary trajectories and 3D genome organization of duplicate genes

This repository contains reproducible scripts for classifying duplicate-gene
expression trajectories and integrating those classifications with 3D genome,
chromatin-accessibility and expression data.

The original trajectory-classification utilities remain in
`evolutionary_trajectory/`. The analysis used for the T2T pear study is
documented in [`study_T2T_pear_3D_duplicate_genes/`](study_T2T_pear_3D_duplicate_genes/README.md).

The T2T study directory includes:

- the four-class two-method consensus analysis;
- the independent CLOUD sensitivity analysis and explicit treatment of the
  CLOUD-specific Sub class;
- hapA 3D-genome integration and stringent hapA-hapB validation;
- source tables for Figure 2, Supplementary Figure S1 and Supplementary Table 17;
- publication figure and manuscript-support scripts.

Run scripts with `-h` where available for command-line options. Paths in the
study workflow are examples and should be adapted to the local installation.
