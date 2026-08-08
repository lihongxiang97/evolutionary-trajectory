#!/usr/bin/env python3
"""Build T2T hapA supplementary figures S1-S7 and S9-S10.

The figures preserve the scientific intent of the original supplement while
using only the final gap-free DS hapA analysis and current trajectory calls.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


CLASS_ORDER = [
    "Conservation",
    "Neofunctionalization(Child)",
    "Neofunctionalization(Parent)",
    "Specialization",
]
CLASS_LABEL = {
    "Conservation": "Con",
    "Neofunctionalization(Child)": "Neo(C)",
    "Neofunctionalization(Parent)": "Neo(P)",
    "Subfunctionalization": "Sub",
    "Specialization": "Spe",
}
CLASS_COLORS = {
    "Conservation": "#4E79A7",
    "Neofunctionalization(Child)": "#F28E2B",
    "Neofunctionalization(Parent)": "#59A14F",
    "Subfunctionalization": "#B07AA1",
    "Specialization": "#E15759",
}
LOOP_COLORS = {
    "Cis-associated": "#4E79A7",
    "Trans-associated": "#E15759",
    "No loop": "#B9B9B9",
}
GGL_COLORS = {
    "GGL-associated": "#59A14F",
    "DGL-associated": "#F28E2B",
    "No loop": "#B9B9B9",
}
REP_COLORS = {"rep1": "#4E79A7", "rep2": "#F28E2B", "merged": "#111111"}


def set_style() -> None:
    sns.set_theme(style="ticks", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.13,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="left",
    )


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "svg", "png", "tiff"):
        kwargs = {
            "bbox_inches": "tight",
            "facecolor": "white",
            "edgecolor": "none",
        }
        if extension in {"png", "tiff"}:
            kwargs["dpi"] = 600
        if extension == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(output_dir / f"{stem}.{extension}", **kwargs)
    plt.close(fig)


def load_contact_decay(source_dir: Path) -> pd.DataFrame:
    tables = []
    for sample in ("rep1", "rep2", "merged"):
        path = source_dir / f"contact_decay.{sample}.tsv"
        frame = pd.read_csv(path, sep="\t")
        frame["sample"] = sample
        frame = frame[
            (frame["distance_bp"] > 0)
            & (frame["mean_contact_per_possible_pair"] > 0)
        ].copy()
        # Remove the sparsely supported extreme-distance tail. At these
        # separations only a handful of chromosome-bin pairs remain and the
        # mean becomes unstable. Retain bins supported by at least 10% of the
        # sample-specific maximum number of possible bin pairs (and >=1,000).
        support_cutoff = max(1000, int(np.ceil(frame["possible_bin_pairs"].max() * 0.10)))
        frame = frame[frame["possible_bin_pairs"] >= support_cutoff].copy()
        frame["support_cutoff_possible_pairs"] = support_cutoff
        tables.append(frame)
    return pd.concat(tables, ignore_index=True)


def make_s1(source_dir: Path, output_dir: Path, source_out: Path) -> dict:
    data = load_contact_decay(source_dir)
    data.to_csv(source_out / "Supplementary_Figure_S1_source.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    for sample in ("rep1", "rep2", "merged"):
        sub = data[data["sample"] == sample].sort_values("distance_bp")
        ax.plot(
            sub["distance_bp"] / 1e6,
            sub["mean_contact_per_possible_pair"],
            color=REP_COLORS[sample],
            linewidth=1.6 if sample == "merged" else 1.1,
            alpha=1.0 if sample == "merged" else 0.82,
            label=sample,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Genomic separation (Mb)")
    ax.set_ylabel("Mean contact frequency")
    ax.legend(frameon=False, title=None)
    sns.despine(ax=ax)
    panel_label(ax, "A")
    save_figure(fig, output_dir, "Supplementary_Figure_S1_contact_decay")

    replicate_rows = []
    pivot = data.pivot_table(
        index="distance_bp",
        columns="sample",
        values="mean_contact_per_possible_pair",
        aggfunc="first",
    ).dropna()
    for sample in ("rep1", "rep2"):
        rho = pivot[[sample, "merged"]].corr(method="spearman").iloc[0, 1]
        replicate_rows.append({"comparison": f"{sample}_vs_merged", "spearman_rho": rho})
    pd.DataFrame(replicate_rows).to_csv(
        source_out / "Supplementary_Figure_S1_summary.tsv", sep="\t", index=False
    )
    return {"S1_rows": int(len(data)), "S1_rep_correlations": replicate_rows}


def load_domain_sizes(source_dir: Path) -> pd.DataFrame:
    rows = []
    for sample in ("rep1", "rep2", "merged"):
        path = source_dir / f"{sample}.primary_domains.bed"
        frame = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            names=["chrom", "start", "end"],
        )
        frame["size_kb"] = (frame["end"] - frame["start"]) / 1000.0
        frame["sample"] = sample
        rows.append(frame)
    return pd.concat(rows, ignore_index=True)


def make_s2(source_dir: Path, output_dir: Path, source_out: Path) -> dict:
    data = load_domain_sizes(source_dir)
    data.to_csv(source_out / "Supplementary_Figure_S2_source.tsv", sep="\t", index=False)

    upper = float(np.nanpercentile(data["size_kb"], 99.5))
    bins = np.linspace(0, upper, 55)
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    for sample in ("rep1", "rep2", "merged"):
        values = data.loc[data["sample"] == sample, "size_kb"].to_numpy()
        hist, edges = np.histogram(values[values <= upper], bins=bins, density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        ax.plot(
            centers,
            hist,
            color=REP_COLORS[sample],
            linewidth=1.6 if sample == "merged" else 1.1,
            label=f"{sample} (n={len(values):,})",
        )
    ax.set_xlabel("TAD size (kb)")
    ax.set_ylabel("Density")
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    panel_label(ax, "A")
    save_figure(fig, output_dir, "Supplementary_Figure_S2_TAD_size")

    summary = (
        data.groupby("sample", sort=False)["size_kb"]
        .agg(n="size", median_kb="median", mean_kb="mean")
        .reset_index()
    )
    summary.to_csv(
        source_out / "Supplementary_Figure_S2_summary.tsv", sep="\t", index=False
    )
    return {"S2_summary": summary.to_dict(orient="records")}


def read_deeptools_profile(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8").splitlines()
    records = []
    for line in lines[2:]:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue
        sample = fields[0]
        values = [float(value) for value in fields[2:] if value != ""]
        positions = np.linspace(-100, 100, len(values))
        records.extend(
            {
                "sample": sample,
                "position_kb": position,
                "ATAC_signal": value,
            }
            for position, value in zip(positions, values)
        )
    return pd.DataFrame(records)


def make_s3(source_dir: Path, output_dir: Path, source_out: Path) -> dict:
    data = read_deeptools_profile(source_dir / "TAD_boundary_ATAC.profile.tsv")
    data.to_csv(source_out / "Supplementary_Figure_S3_source.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    for sample in ("rep1", "rep2"):
        sub = data[data["sample"] == sample]
        ax.plot(
            sub["position_kb"],
            sub["ATAC_signal"],
            color=REP_COLORS[sample],
            linewidth=1.5,
            label=sample,
        )
    ax.axvline(0, color="#666666", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Distance from TAD boundary (kb)")
    ax.set_ylabel("ATAC signal (RPGC)")
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    panel_label(ax, "A")
    save_figure(fig, output_dir, "Supplementary_Figure_S3_TAD_boundary_ATAC")

    central = data[data["position_kb"].abs() <= 10].groupby("sample")[
        "ATAC_signal"
    ].mean()
    flank = data[data["position_kb"].abs() >= 70].groupby("sample")[
        "ATAC_signal"
    ].mean()
    summary = pd.DataFrame(
        {
            "sample": central.index,
            "boundary_mean": central.values,
            "distal_flank_mean": flank.reindex(central.index).values,
        }
    )
    summary["boundary_to_flank_ratio"] = (
        summary["boundary_mean"] / summary["distal_flank_mean"]
    )
    summary.to_csv(
        source_out / "Supplementary_Figure_S3_summary.tsv", sep="\t", index=False
    )
    return {"S3_summary": summary.to_dict(orient="records")}


def violin_box(
    ax: plt.Axes,
    data: pd.DataFrame,
    x: str,
    y: str,
    order: list[str],
    palette: dict[str, str],
) -> None:
    sns.violinplot(
        data=data,
        x=x,
        y=y,
        order=order,
        palette=palette,
        inner=None,
        linewidth=0.7,
        cut=0,
        density_norm="width",
        ax=ax,
    )
    sns.boxplot(
        data=data,
        x=x,
        y=y,
        order=order,
        width=0.22,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.7},
        medianprops={"color": "#111111", "linewidth": 1.0},
        whiskerprops={"color": "#333333", "linewidth": 0.7},
        capprops={"color": "#333333", "linewidth": 0.7},
        ax=ax,
    )


def make_s4_s5(
    expression_path: Path,
    loop_dir: Path,
    output_dir: Path,
    source_out: Path,
) -> dict:
    expression = pd.read_csv(expression_path, sep="\t")
    expression["log2_TPM_plus_1"] = np.log2(expression["TPM"].clip(lower=0) + 1)
    cis_genes = set(
        pd.read_csv(loop_dir / "cis.genes.txt", header=None)[0].astype(str)
    )
    trans_genes = set(
        pd.read_csv(loop_dir / "trans.genes.txt", header=None)[0].astype(str)
    )

    s4_frames = []
    for label, mask in (
        ("Cis-associated", expression["gene"].isin(cis_genes)),
        ("Trans-associated", expression["gene"].isin(trans_genes)),
        ("No loop", ~expression["any_loop"].astype(bool)),
    ):
        part = expression.loc[mask, ["gene", "TPM", "log2_TPM_plus_1"]].copy()
        part["group"] = label
        s4_frames.append(part)
    s4 = pd.concat(s4_frames, ignore_index=True)
    s4.to_csv(source_out / "Supplementary_Figure_S4_source.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    s4_order = ["Cis-associated", "Trans-associated", "No loop"]
    violin_box(ax, s4, "group", "log2_TPM_plus_1", s4_order, LOOP_COLORS)
    ax.set_xlabel("")
    ax.set_ylabel(r"$\log_2$(fruit TPM + 1)")
    ax.set_xticklabels(
        [
            f"Cis-associated\n(n={(s4['group'] == 'Cis-associated').sum():,})",
            f"Trans-associated\n(n={(s4['group'] == 'Trans-associated').sum():,})",
            f"No loop\n(n={(s4['group'] == 'No loop').sum():,})",
        ]
    )
    sns.despine(ax=ax)
    panel_label(ax, "A")
    save_figure(fig, output_dir, "Supplementary_Figure_S4_loop_expression")

    s5_frames = []
    for label, mask in (
        ("GGL-associated", expression["GGL"].astype(bool)),
        ("DGL-associated", expression["DGL"].astype(bool)),
        ("No loop", ~expression["any_loop"].astype(bool)),
    ):
        part = expression.loc[mask, ["gene", "TPM", "log2_TPM_plus_1"]].copy()
        part["group"] = label
        s5_frames.append(part)
    s5 = pd.concat(s5_frames, ignore_index=True)
    s5.to_csv(source_out / "Supplementary_Figure_S5_source.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    s5_order = ["GGL-associated", "DGL-associated", "No loop"]
    violin_box(ax, s5, "group", "log2_TPM_plus_1", s5_order, GGL_COLORS)
    ax.set_xlabel("")
    ax.set_ylabel(r"$\log_2$(fruit TPM + 1)")
    ax.set_xticklabels(
        [
            f"GGL-associated\n(n={(s5['group'] == 'GGL-associated').sum():,})",
            f"DGL-associated\n(n={(s5['group'] == 'DGL-associated').sum():,})",
            f"No loop\n(n={(s5['group'] == 'No loop').sum():,})",
        ]
    )
    sns.despine(ax=ax)
    panel_label(ax, "A")
    save_figure(fig, output_dir, "Supplementary_Figure_S5_GGL_DGL_expression")

    summary_rows = []
    for figure_id, data in (("S4", s4), ("S5", s5)):
        for group, sub in data.groupby("group", sort=False):
            summary_rows.append(
                {
                    "figure": figure_id,
                    "group": group,
                    "n": len(sub),
                    "median_TPM": float(sub["TPM"].median()),
                    "mean_log2_TPM_plus_1": float(sub["log2_TPM_plus_1"].mean()),
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        source_out / "Supplementary_Figures_S4_S5_summary.tsv",
        sep="\t",
        index=False,
    )
    return {"S4_S5_summary": summary.to_dict(orient="records")}


def make_s6(cloud_path: Path, output_dir: Path, source_out: Path) -> dict:
    data = pd.read_csv(cloud_path, sep="\t")
    probability_column = {
        "Conservation": "prob_Conservation",
        "Neofunctionalization(Parent)": "prob_Neofunctionalization_Parent",
        "Neofunctionalization(Child)": "prob_Neofunctionalization_Child",
        "Subfunctionalization": "prob_Subfunctionalization",
        "Specialization": "prob_Specialization",
    }
    data["assigned_class_probability"] = [
        row[probability_column[row["CLOUD_class"]]]
        for _, row in data.iterrows()
    ]
    cloud_order = [
        "Conservation",
        "Neofunctionalization(Child)",
        "Neofunctionalization(Parent)",
        "Subfunctionalization",
        "Specialization",
    ]
    data["class_label"] = data["CLOUD_class"].map(CLASS_LABEL)
    data.to_csv(source_out / "Supplementary_Figure_S6_source.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(5.7, 3.7))
    violin_box(
        ax,
        data,
        "CLOUD_class",
        "assigned_class_probability",
        cloud_order,
        CLASS_COLORS,
    )
    ax.set_ylim(0, 1.03)
    ax.set_xlabel("")
    ax.set_ylabel("Assigned-class posterior probability")
    ax.set_xticklabels(
        [
            f"{CLASS_LABEL[item]}\n(n={(data['CLOUD_class'] == item).sum():,})"
            for item in cloud_order
        ]
    )
    sns.despine(ax=ax)
    panel_label(ax, "A")
    save_figure(fig, output_dir, "Supplementary_Figure_S6_CLOUD_confidence")

    summary = (
        data.groupby("CLOUD_class", sort=False)["assigned_class_probability"]
        .agg(n="size", median="median", mean="mean")
        .reset_index()
    )
    summary.to_csv(
        source_out / "Supplementary_Figure_S6_summary.tsv", sep="\t", index=False
    )
    return {"S6_summary": summary.to_dict(orient="records")}


def normalize_confusion(
    data: pd.DataFrame, row_col: str, col_col: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    allowed = set(CLASS_ORDER)
    sub = data[data[row_col].isin(allowed) & data[col_col].isin(allowed)].copy()
    counts = pd.crosstab(sub[row_col], sub[col_col]).reindex(
        index=CLASS_ORDER, columns=CLASS_ORDER, fill_value=0
    )
    fractions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)
    return counts, fractions


def make_s7(
    trajectory_path: Path,
    cloud_path: Path,
    output_dir: Path,
    source_out: Path,
) -> dict:
    trajectory = pd.read_csv(trajectory_path, sep="\t")
    cloud = pd.read_csv(cloud_path, sep="\t")[
        ["duplication_type", "parent", "child", "CLOUD_class"]
    ]
    data = trajectory.merge(
        cloud,
        on=["duplication_type", "parent", "child"],
        how="inner",
        validate="one_to_one",
    )
    comparisons = [
        ("binary_TPM1_class", "CDROM_class", "Localization vs CDROM"),
        ("binary_TPM1_class", "CLOUD_class", "Localization vs CLOUD"),
        ("CDROM_class", "CLOUD_class", "CDROM vs CLOUD"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.45), constrained_layout=True)
    source_frames = []
    summaries = []
    for panel, (ax, (rows, cols, title)) in enumerate(
        zip(axes, comparisons), start=1
    ):
        counts, fractions = normalize_confusion(data, rows, cols)
        annotations = counts.astype(str) + "\n" + fractions.applymap(
            lambda value: "" if pd.isna(value) else f"{100 * value:.0f}%"
        )
        sns.heatmap(
            fractions,
            vmin=0,
            vmax=1,
            cmap="Blues",
            square=True,
            linewidths=0.4,
            linecolor="white",
            annot=annotations,
            fmt="",
            cbar=panel == 3,
            cbar_kws={"label": "Row fraction", "shrink": 0.72},
            ax=ax,
        )
        ax.set_title(title, pad=7)
        ax.set_xlabel(CLASS_LABEL.get(cols, cols))
        ax.set_ylabel(CLASS_LABEL.get(rows, rows))
        ax.set_xticklabels([CLASS_LABEL[item] for item in CLASS_ORDER], rotation=45)
        ax.set_yticklabels([CLASS_LABEL[item] for item in CLASS_ORDER], rotation=0)
        panel_label(ax, chr(64 + panel))
        long_counts = counts.stack().rename("n").reset_index()
        long_fractions = fractions.stack(dropna=False).rename("row_fraction").reset_index()
        long = long_counts.merge(long_fractions, on=[rows, cols], how="outer")
        long["comparison"] = title
        source_frames.append(long)
        agreement = (
            np.trace(counts.to_numpy()) / counts.to_numpy().sum()
            if counts.to_numpy().sum()
            else np.nan
        )
        summaries.append(
            {
                "comparison": title,
                "n": int(counts.to_numpy().sum()),
                "exact_agreement": float(agreement),
            }
        )
    save_figure(fig, output_dir, "Supplementary_Figure_S7_method_confusion")
    pd.concat(source_frames, ignore_index=True).to_csv(
        source_out / "Supplementary_Figure_S7_source.tsv", sep="\t", index=False
    )
    pd.DataFrame(summaries).to_csv(
        source_out / "Supplementary_Figure_S7_summary.tsv", sep="\t", index=False
    )
    return {"S7_summary": summaries}


def make_s9_s10(
    tau_path: Path,
    output_dir: Path,
    source_out: Path,
) -> dict:
    data = pd.read_csv(tau_path, sep="\t")
    data = data[data["two_method_agreement"].isin(CLASS_ORDER)].copy()
    data.to_csv(source_out / "Supplementary_Figure_S9_source.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(
        1, 3, figsize=(10.1, 3.7), sharey=True, constrained_layout=True
    )
    for panel, (ax, duplication_type) in enumerate(
        zip(axes, ["TD", "PD", "TRD"]), start=1
    ):
        sub = data[data["duplication_type"] == duplication_type]
        violin_box(
            ax,
            sub,
            "two_method_agreement",
            "ancestor_tau",
            CLASS_ORDER,
            CLASS_COLORS,
        )
        ax.set_ylim(-0.03, 1.04)
        ax.set_xlabel("")
        ax.set_ylabel("Ancestral ortholog tissue specificity (τ)")
        ax.set_title(duplication_type)
        ax.set_xticklabels(
            [
                f"{CLASS_LABEL[item]}\n(n={(sub['two_method_agreement'] == item).sum():,})"
                for item in CLASS_ORDER
            ],
            rotation=0,
        )
        sns.despine(ax=ax)
        panel_label(ax, chr(64 + panel))
    save_figure(fig, output_dir, "Supplementary_Figure_S9_ancestral_tau")

    s9_summary = (
        data.groupby(["duplication_type", "two_method_agreement"], sort=False)[
            "ancestor_tau"
        ]
        .agg(n="size", median="median", mean="mean")
        .reset_index()
    )
    s9_summary.to_csv(
        source_out / "Supplementary_Figure_S9_summary.tsv", sep="\t", index=False
    )

    data.to_csv(source_out / "Supplementary_Figure_S10_source.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 7.0), constrained_layout=True)
    rng = np.random.default_rng(20260729)
    s10_summary = []
    for panel, (ax, trajectory) in enumerate(zip(axes.flat, CLASS_ORDER), start=1):
        sub = data[data["two_method_agreement"] == trajectory].copy()
        if len(sub) > 2500:
            selected = rng.choice(sub.index.to_numpy(), size=2500, replace=False)
            plot_sub = sub.loc[selected]
        else:
            plot_sub = sub
        ax.scatter(
            plot_sub["ancestor_tau"],
            plot_sub["parent_tau"],
            s=8,
            alpha=0.25,
            linewidths=0,
            color="#4E79A7",
            label="Parent",
            rasterized=True,
        )
        ax.scatter(
            plot_sub["ancestor_tau"],
            plot_sub["child_tau"],
            s=8,
            alpha=0.25,
            linewidths=0,
            color="#E15759",
            label="Child",
            rasterized=True,
        )
        ax.plot([0, 1], [0, 1], color="#555555", linewidth=0.8, linestyle="--")
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{CLASS_LABEL[trajectory]} (n={len(sub):,})")
        ax.set_xlabel("Ancestral ortholog τ")
        ax.set_ylabel("Duplicate-gene τ")
        if panel == 1:
            ax.legend(frameon=False, loc="lower right")
        sns.despine(ax=ax)
        panel_label(ax, chr(64 + panel))
        for duplicate_member, column in (("Parent", "parent_tau"), ("Child", "child_tau")):
            rho = sub[["ancestor_tau", column]].corr(method="spearman").iloc[0, 1]
            s10_summary.append(
                {
                    "trajectory": trajectory,
                    "duplicate_member": duplicate_member,
                    "n": len(sub),
                    "spearman_rho": float(rho),
                }
            )
    save_figure(fig, output_dir, "Supplementary_Figure_S10_duplicate_ancestor_tau")
    pd.DataFrame(s10_summary).to_csv(
        source_out / "Supplementary_Figure_S10_summary.tsv", sep="\t", index=False
    )
    return {
        "S9_summary": s9_summary.to_dict(orient="records"),
        "S10_summary": s10_summary,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hapA-source", required=True, type=Path)
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--cloud", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--tau", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--source-output-dir", required=True, type=Path)
    args = parser.parse_args()

    set_style()
    args.source_output_dir.mkdir(parents=True, exist_ok=True)
    summaries = {}
    summaries.update(
        make_s1(args.hapA_source, args.output_dir, args.source_output_dir)
    )
    summaries.update(
        make_s2(args.hapA_source, args.output_dir, args.source_output_dir)
    )
    summaries.update(
        make_s3(args.hapA_source, args.output_dir, args.source_output_dir)
    )
    summaries.update(
        make_s4_s5(
            args.expression,
            args.hapA_source / "gene_loops",
            args.output_dir,
            args.source_output_dir,
        )
    )
    summaries.update(
        make_s6(args.cloud, args.output_dir, args.source_output_dir)
    )
    summaries.update(
        make_s7(
            args.trajectory,
            args.cloud,
            args.output_dir,
            args.source_output_dir,
        )
    )
    summaries.update(
        make_s9_s10(args.tau, args.output_dir, args.source_output_dir)
    )
    (args.source_output_dir / "supplementary_figures_S1_S10_summary.json").write_text(
        json.dumps(summaries, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
