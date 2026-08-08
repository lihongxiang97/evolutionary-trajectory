#!/usr/bin/env python3
from pathlib import Path
from bisect import bisect_right
import json
import numpy as np
import pandas as pd
import pyBigWig

EVOL = Path(__file__).resolve().parents[1]
PROJECT = EVOL.parent
GENOME3D = PROJECT / "05.3D_genome"
OUT = Path(__file__).resolve().parent / "source"
OUT.mkdir(parents=True, exist_ok=True)

TRAJ = EVOL / "08.trajectories/results/trajectory_classifications.tsv"
SINGLE = EVOL / "06.kaks/Pbr_Ppe_singlecopy/pairs.tsv"
TSS = GENOME3D / "00.metadata/hapA.TSS.bed"
GENE3D = GENOME3D / "07.integration/gene_3D_expression.tsv"
TADS = GENOME3D / "03.TAD/merged.primary_domains.bed"
ATAC = [
    PROJECT / "02.ATAC-seq/bigwig/ds15-ATAC-rep1.RPGC.bw",
    PROJECT / "02.ATAC-seq/bigwig/ds15-ATAC-rep2.RPGC.bw",
]

CLASSES = [
    "Conservation",
    "Neofunctionalization(Child)",
    "Neofunctionalization(Parent)",
    "Specialization",
]
SHORT = {
    "Conservation": "Con",
    "Neofunctionalization(Child)": "Neo(C)",
    "Neofunctionalization(Parent)": "Neo(P)",
    "Specialization": "Spe",
}

traj = pd.read_csv(TRAJ, sep="\t")
traj = traj[
    (traj["analysis_status"] == "classified")
    & traj["two_method_agreement"].isin(CLASSES)
].copy()
traj["trajectory"] = traj["two_method_agreement"].map(SHORT)
traj["pair_id"] = (
    traj["duplication_type"] + "|" + traj["parent"] + "|" + traj["child"]
)

records = []
for row in traj.itertuples(index=False):
    records.append(
        {
            "trajectory": row.trajectory,
            "duplication_type": row.duplication_type,
            "role": "P",
            "gene": row.parent,
            "pair_id": row.pair_id,
            "reference_definition": "duplicate_parent",
        }
    )
    records.append(
        {
            "trajectory": row.trajectory,
            "duplication_type": row.duplication_type,
            "role": "C",
            "gene": row.child,
            "pair_id": row.pair_id,
            "reference_definition": "duplicate_child",
        }
    )

single = pd.read_csv(SINGLE, sep="\t", header=None, names=["pear_gene", "peach_gene"])
single_genes = sorted(single["pear_gene"].dropna().astype(str).unique())
# Single is a pear chromatin reference, not the peach ancestor itself:
# pear genes retained as one-to-one single-copy orthologs with peach.
for trajectory in map(SHORT.get, CLASSES):
    for gene in single_genes:
        records.append(
            {
                "trajectory": trajectory,
                "duplication_type": "Reference",
                "role": "S",
                "gene": gene,
                "pair_id": "",
                "reference_definition": "pear_peach_singlecopy_reference",
            }
        )
assign = pd.DataFrame(records)

tss = pd.read_csv(
    TSS,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "gene", "score", "strand"],
)
tss["tss"] = tss["start"].astype(int)
tss = tss.drop_duplicates("gene")
tss_map = tss.set_index("gene")[["chrom", "tss", "strand"]]

gene3d = pd.read_csv(GENE3D, sep="\t")
for column in ["any_loop", "GGL", "DGL"]:
    gene3d[column] = gene3d[column].astype(str).str.lower().eq("true")
gene3d = gene3d.drop_duplicates("gene").set_index("gene")

domains = pd.read_csv(
    TADS,
    sep="\t",
    header=None,
    usecols=[0, 1, 2],
    names=["chrom", "start", "end"],
)
domain_by_chrom = {}
for chrom, block in domains.groupby("chrom"):
    intervals = sorted(zip(block["start"].astype(int), block["end"].astype(int)))
    domain_by_chrom[chrom] = (intervals, [start for start, _ in intervals])


def tad_location(chrom, position, boundary_window=20000):
    if chrom not in domain_by_chrom:
        return "Outside"
    intervals, starts = domain_by_chrom[chrom]
    idx = bisect_right(starts, position)
    candidates = intervals[max(0, idx - 2) : min(len(intervals), idx + 2)]
    for start, end in candidates:
        if abs(position - start) <= boundary_window or abs(position - end) <= boundary_window:
            return "Boundary"
    for start, end in candidates:
        if start < position < end:
            return "Interior"
    return "Outside"


unique_genes = sorted(set(assign["gene"]) & set(tss_map.index))
bws = [pyBigWig.open(str(path)) for path in ATAC]
chrom_sizes = bws[0].chroms()
window = 3000
n_bins = 60
bin_centers = np.linspace(-window + window / n_bins, window - window / n_bins, n_bins)
profile_rows = []
feature_rows = []

for gene in unique_genes:
    chrom, position, strand = tss_map.loc[gene]
    position = int(position)
    values = []
    if chrom in chrom_sizes:
        start = max(0, position - window)
        end = min(chrom_sizes[chrom], position + window)
        if end > start:
            for bw in bws:
                signal = bw.stats(chrom, start, end, nBins=n_bins, exact=True)
                signal = np.array(
                    [np.nan if value is None else value for value in signal],
                    dtype=float,
                )
                if strand == "-":
                    signal = signal[::-1]
                values.append(signal)
    if values:
        profile = np.nanmean(np.vstack(values), axis=0)
    else:
        profile = np.full(n_bins, np.nan)
    mean_tss_atac = float(np.nanmean(profile[20:40])) if np.isfinite(profile).any() else np.nan
    for center, value in zip(bin_centers, profile):
        profile_rows.append({"gene": gene, "position": center, "ATAC_RPGC": value})

    info = gene3d.loc[gene] if gene in gene3d.index else None
    feature_rows.append(
        {
            "gene": gene,
            "chrom": chrom,
            "tss": position,
            "strand": strand,
            "mean_ATAC_TSS_1kb": mean_tss_atac,
            "PC1": info["PC1"] if info is not None else np.nan,
            "compartment": info["compartment"] if info is not None else "Unknown",
            "TAD_location": tad_location(chrom, position),
            "any_loop": bool(info["any_loop"]) if info is not None else False,
            "GGL": bool(info["GGL"]) if info is not None else False,
            "DGL": bool(info["DGL"]) if info is not None else False,
        }
    )

for bw in bws:
    bw.close()

profiles = pd.DataFrame(profile_rows)
features = pd.DataFrame(feature_rows)
annotated = assign.merge(features, on="gene", how="left")
annotated.to_csv(OUT / "Figure4_gene_features.tsv.gz", sep="\t", index=False)

profile_assign = assign[["trajectory", "role", "gene"]].drop_duplicates()
profile_long = profile_assign.merge(profiles, on="gene", how="inner")
summary = (
    profile_long.groupby(["trajectory", "role", "position"])["ATAC_RPGC"]
    .agg(["mean", "std", "count"])
    .reset_index()
)
summary["sem"] = summary["std"] / np.sqrt(summary["count"])
summary["ci95_low"] = summary["mean"] - 1.96 * summary["sem"]
summary["ci95_high"] = summary["mean"] + 1.96 * summary["sem"]
summary.to_csv(OUT / "Figure4_ATAC_profile_summary.tsv.gz", sep="\t", index=False)

qc = {
    "trajectory_pairs": int(len(traj)),
    "singlecopy_reference_genes": int(len(single_genes)),
    "unique_genes_with_tss": int(len(unique_genes)),
    "role_records": assign.groupby(["trajectory", "role"]).size().to_dict(),
    "reference_note": (
        "Single is a pear-peach one-to-one single-copy pear-gene chromatin reference; "
        "no peach Hi-C/ATAC data are used."
    ),
    "TAD_boundary_window_bp": 20000,
    "ATAC_window_bp": 3000,
    "ATAC_bin_bp": 100,
}
qc["role_records"] = {"|".join(key): value for key, value in qc["role_records"].items()}
(OUT / "Figure4_preparation_summary.json").write_text(
    json.dumps(qc, indent=2), encoding="utf-8"
)
print(json.dumps(qc, indent=2))



