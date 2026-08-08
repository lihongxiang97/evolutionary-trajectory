from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


COLORS = {
    "TD": "#0072B2",
    "PD": "#D55E00",
    "TRD": "#009E73",
    "Con": "#4C78A8",
    "Neo(C)": "#E45756",
    "Neo(P)": "#F2CF5B",
    "Spe": "#72B7B2",
    "Young": "#D55E00",
    "Old": "#4C78A8",
}
TRAJ = ["Con", "Neo(C)", "Neo(P)", "Spe"]
TYPES = ["TD", "PD", "TRD"]
ROLES = [("Con", "P"), ("Con", "C"), ("Neo(C)", "P"), ("Neo(C)", "C"),
         ("Neo(P)", "P"), ("Neo(P)", "C"), ("Spe", "P"), ("Spe", "C")]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.5,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.15, y: float = 1.08) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left")


def clean_ax(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def export(fig: plt.Figure, outdir: Path, stem: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def fig_s13(indir: Path, outdir: Path) -> None:
    comp = pd.read_csv(indir / "consensus_composition_with_CI.tsv", sep="\t")
    ret = pd.read_csv(indir / "consensus_retention.tsv", sep="\t")
    assoc = pd.read_csv(indir / "consensus_duplication_type_association.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 3.45))
    gs = fig.add_gridspec(1, 3, width_ratios=[2.1, 1.2, 1.25], wspace=0.55)

    ax = fig.add_subplot(gs[0, 0])
    positions, labels = [], []
    pos = 0
    for aset, method_label in [("two_method_primary", "Two-method"),
                               ("three_method", "Three-method")]:
        for dtype in TYPES:
            subset = comp[(comp.analysis_set == aset) & (comp.duplication_type == dtype)]
            bottom = 0
            for tr in TRAJ:
                row = subset[subset.trajectory == tr]
                value = float(row.percent.iloc[0]) if len(row) else 0
                ax.bar(pos, value, bottom=bottom, width=0.72, color=COLORS[tr],
                       edgecolor="white", linewidth=0.35)
                bottom += value
            positions.append(pos)
            labels.append(f"{dtype}\n{method_label}")
            pos += 1
        pos += 0.55
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylim(0, 100)
    ax.set_ylabel("Trajectory composition (%)")
    clean_ax(ax)
    fig.legend(handles=[Patch(facecolor=COLORS[t], label=t) for t in TRAJ],
               frameon=False, ncol=4, loc="upper left",
               bbox_to_anchor=(0.08, 0.985), columnspacing=1.0)
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    r = ret.set_index("duplication_type").loc[["All", *TYPES]].reset_index()
    x = np.arange(len(r))
    y = r.strict_retained_percent.to_numpy()
    lo = y - r.ci95_low_percent.to_numpy()
    hi = r.ci95_high_percent.to_numpy() - y
    colors = ["#666666", *[COLORS[t] for t in TYPES]]
    ax.bar(x, y, width=0.65, color=colors, edgecolor="none")
    ax.errorbar(x, y, yerr=[lo, hi], fmt="none", ecolor="black",
                capsize=2, lw=0.8)
    for i, row in r.iterrows():
        ax.text(i, y[i] + hi[i] + 2, f"{int(row.strict_retained_n)}/{int(row.primary_n)}",
                ha="center", va="bottom", fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels(r.duplication_type)
    ax.set_ylim(0, max(35, (y + hi).max() + 9))
    ax.set_ylabel("Retained by CLOUD (%)")
    clean_ax(ax)
    panel_label(ax, "B")

    ax = fig.add_subplot(gs[0, 2])
    a = assoc.copy()
    y = np.arange(len(a))[::-1]
    for yi, (_, row) in zip(y, a.iterrows()):
        color = "#333333" if row.analysis_set == "two_method_primary" else "#7A5195"
        ax.errorbar(row.cramers_V, yi,
                    xerr=[[row.cramers_V - row.ci95_low],
                          [row.ci95_high - row.cramers_V]],
                    fmt="o", color=color, capsize=2, lw=1, ms=4)
    ax.set_yticks(y)
    ax.set_yticklabels(["Two-method primary", "Three-method sensitivity"])
    ax.set_xlabel("Cramér's V (95% CI)")
    ax.set_xlim(0, max(0.5, assoc.ci95_high.max() + 0.06))
    ax.axvline(0, color="#999999", lw=0.7)
    clean_ax(ax)
    panel_label(ax, "C", x=-0.25)

    fig.subplots_adjust(left=0.10, right=0.98, top=0.82, bottom=0.24)
    export(fig, outdir, "Supplementary_Figure_S13_consensus_robustness")


def forest_young_old(ax: plt.Axes, data: pd.DataFrame, outcome: str,
                     ylabel: str, odds: bool = False) -> None:
    d = data[data.outcome == outcome].copy()
    offsets = {"all_annotated": -0.12, "complete_ORF_expressed": 0.12}
    markers = {"all_annotated": "o", "complete_ORF_expressed": "s"}
    for i, dtype in enumerate(TYPES):
        for spec in offsets:
            row = d[(d.duplication_type == dtype) & (d.sensitivity == spec)].iloc[0]
            x = row.effect_young_vs_old
            lo, hi = row.ci95_low, row.ci95_high
            ax.errorbar(x, i + offsets[spec], xerr=[[x - lo], [hi - x]],
                        fmt=markers[spec], color=COLORS[dtype], ms=4,
                        capsize=2, lw=0.9, markeredgecolor="white",
                        markeredgewidth=0.35)
    ax.set_yticks(range(3))
    ax.set_yticklabels(TYPES)
    ax.set_xlabel(ylabel)
    ax.axvline(1 if odds else 0, color="#777777", lw=0.75, ls="--")
    if odds:
        ax.set_xscale("log")
    clean_ax(ax)


def fig_s14(indir: Path, outdir: Path) -> None:
    d = pd.read_csv(indir / "adjusted_young_old_chromatin_models.tsv", sep="\t")
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0))
    specs = [
        ("log2_TPM_plus_1", "Adjusted β: log$_2$(TPM + 1)", False),
        ("mean_ATAC_TSS_2kb", "Adjusted β: ATAC signal", False),
        ("PC1", "Adjusted β: compartment PC1", False),
        ("log2_contact_CPM_plus_1", "Adjusted β: log$_2$(contact CPM + 1)", False),
        ("compartment_A", "Adjusted odds ratio: A compartment", True),
        ("TAD_boundary", "Adjusted odds ratio: TAD boundary", True),
    ]
    for idx, (ax, (outcome, label, odds)) in enumerate(zip(axes.flat, specs)):
        forest_young_old(ax, d, outcome, label, odds)
        panel_label(ax, chr(ord("A") + idx), x=-0.2)
    handles = [
        Line2D([0], [0], marker="o", color="#555555", lw=0, label="All annotated",
               markerfacecolor="#555555", markersize=4),
        Line2D([0], [0], marker="s", color="#555555", lw=0,
               label="Complete ORF, expressed, conflict-excluded",
               markerfacecolor="#555555", markersize=4),
    ]
    fig.legend(handles=handles, frameon=False, loc="lower center",
               bbox_to_anchor=(0.5, -0.015), ncol=2)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.98, bottom=0.12,
                        wspace=0.55, hspace=0.48)
    export(fig, outdir, "Supplementary_Figure_S14_adjusted_age_effects")


def fig_s15(indir: Path, outdir: Path, te_path: Path | None = None) -> None:
    q = pd.read_csv(indir / "young_old_annotation_expression_QC.tsv", sep="\t")
    cont_features = [
        ("gene_length", "Gene length"),
        ("CDS_length", "CDS length"),
        ("intron_count", "Intron count"),
        ("six_tissue_mean_TPM", "Six-tissue mean TPM"),
    ]
    binary_features = [
        ("intronless", "Intronless"),
        ("expressed_any_TPM1", "Expressed in ≥1 tissue"),
        ("fruit_expressed_TPM1", "Expressed in fruit"),
    ]

    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(2, 2, wspace=0.5, hspace=0.72)

    ax = fig.add_subplot(gs[0, 0])
    ybase = np.arange(len(cont_features))[::-1]
    for j, dtype in enumerate(TYPES):
        sub = q[q.duplication_type == dtype].set_index("feature")
        for yi, (feature, _) in zip(ybase, cont_features):
            row = sub.loc[feature]
            ax.plot(row.effect, yi + (j - 1) * 0.16, "o", color=COLORS[dtype],
                    ms=4)
    ax.axvline(0, color="#777777", lw=0.75, ls="--")
    ax.set_yticks(ybase)
    ax.set_yticklabels([label for _, label in cont_features])
    ax.set_xlabel("Rank-biserial effect (Young vs Old)")
    clean_ax(ax)
    panel_label(ax, "A", x=-0.22)

    ax = fig.add_subplot(gs[0, 1])
    ybase = np.arange(len(binary_features))[::-1]
    for j, dtype in enumerate(TYPES):
        sub = q[q.duplication_type == dtype].set_index("feature")
        for yi, (feature, _) in zip(ybase, binary_features):
            row = sub.loc[feature]
            ax.plot(row.effect, yi + (j - 1) * 0.16, "o", color=COLORS[dtype],
                    ms=4)
    ax.axvline(1, color="#777777", lw=0.75, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(ybase)
    ax.set_yticklabels([label for _, label in binary_features])
    ax.set_xlabel("Odds ratio (Young vs Old)")
    clean_ax(ax)
    panel_label(ax, "B", x=-0.22)

    ax = fig.add_subplot(gs[1, 0])
    features = ["annotation_complete", "expressed_any_TPM1", "fruit_expressed_TPM1"]
    labels = ["Complete annotation", "Expressed in ≥1 tissue", "Expressed in fruit"]
    trd = q[q.duplication_type == "TRD"].set_index("feature")
    x = np.arange(len(features))
    w = 0.36
    young = [100 * trd.loc[f, "young_summary"] for f in features]
    old = [100 * trd.loc[f, "old_summary"] for f in features]
    ax.bar(x - w / 2, young, width=w, color=COLORS["Young"], label="Young")
    ax.bar(x + w / 2, old, width=w, color=COLORS["Old"], label="Old")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("TRD genes (%)")
    ax.set_ylim(0, 108)
    clean_ax(ax)
    panel_label(ax, "C", x=-0.22)

    ax = fig.add_subplot(gs[1, 1])
    if te_path and te_path.exists():
        te = pd.read_csv(te_path, sep="\t")
        trd_te = te[te.duplication_type == "TRD"].copy()
        categories = ["gene_body_TE", "promoter_TE"]
        labels = ["Gene-body TE overlap", "Promoter TE overlap"]
        x = np.arange(len(categories))
        w = 0.36
        yv, ov = [], []
        for feature in categories:
            row = trd_te[trd_te.feature == feature].iloc[0]
            yv.append(100 * row.young_fraction)
            ov.append(100 * row.old_fraction)
        ax.bar(x - w / 2, yv, width=w, color=COLORS["Young"], label="Young")
        ax.bar(x + w / 2, ov, width=w, color=COLORS["Old"], label="Old")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=18, ha="right")
        ax.set_ylabel("TRD genes (%)")
    else:
        ax.text(0.5, 0.55, "TE sensitivity pending", ha="center", va="center",
                transform=ax.transAxes, fontsize=9)
        ax.text(0.5, 0.42, "Panel is populated after RepeatMasker completes",
                ha="center", va="center", transform=ax.transAxes, fontsize=7,
                color="#666666")
        ax.set_xticks([])
        ax.set_yticks([])
    clean_ax(ax)
    panel_label(ax, "D", x=-0.22)

    age_handles = [
        Line2D([0], [0], marker="s", color="none",
               markerfacecolor=COLORS[age], markeredgecolor="none",
               label=age, markersize=5)
        for age in ["Young", "Old"]
    ]
    fig.legend(handles=age_handles, frameon=False, ncol=2, loc="center",
               bbox_to_anchor=(0.5, 0.50))

    handles = [Line2D([0], [0], marker="o", color=COLORS[t], lw=0,
                      markerfacecolor=COLORS[t], label=t, markersize=4)
               for t in TYPES]
    fig.legend(handles=handles, frameon=False, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.005))
    fig.subplots_adjust(left=0.14, right=0.98, top=0.98, bottom=0.14)
    export(fig, outdir, "Supplementary_Figure_S15_young_gene_QC")


def fig_s16(indir: Path, outdir: Path) -> None:
    m = pd.read_csv(indir / "matched_single_copy_chromatin_tests.tsv", sep="\t")
    b = pd.read_csv(indir / "matched_single_copy_covariate_balance.tsv", sep="\t")
    loop = pd.read_csv(indir / "adjusted_loop_expression_models.tsv", sep="\t")
    m = m[m.specification == "genomic_plus_expression"].copy()

    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(3, 2, hspace=0.6, wspace=0.52)
    outcome_specs = [
        ("mean_ATAC_TSS_1kb", "Paired median difference in ATAC signal", False),
        ("compartment_A", "Matched odds ratio: A compartment", True),
        ("TAD_boundary", "Matched odds ratio: TAD boundary", True),
        ("any_loop", "Matched odds ratio: loop association", True),
    ]
    ylabels = [f"{t} {r}" for t, r in ROLES][::-1]
    for idx, (outcome, xlabel, odds) in enumerate(outcome_specs):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        d = m[m.outcome == outcome].set_index(["trajectory", "role"])
        for yi, (tr, role) in enumerate(ROLES[::-1]):
            row = d.loc[(tr, role)]
            x, lo, hi = row.effect, row.ci95_low, row.ci95_high
            ax.errorbar(x, yi, xerr=[[x - lo], [hi - x]], fmt="o",
                        color=COLORS[tr], ms=4, capsize=2, lw=0.9)
        ax.axvline(1 if odds else 0, color="#777777", lw=0.75, ls="--")
        if odds:
            ax.set_xscale("log")
        ax.set_yticks(range(len(ylabels)))
        ax.set_yticklabels(ylabels)
        ax.set_xlabel(xlabel)
        clean_ax(ax)
        panel_label(ax, chr(ord("A") + idx), x=-0.23)

    ax = fig.add_subplot(gs[2, 0])
    bb = b[b.specification == "genomic_plus_expression"].copy()
    max_smd = (bb.assign(abs_smd=bb.SMD_after_matching.abs())
                 .groupby(["trajectory", "role"], as_index=False).abs_smd.max())
    max_smd["label"] = max_smd.trajectory + " " + max_smd.role
    max_smd = max_smd.set_index(["trajectory", "role"]).loc[ROLES].reset_index()
    ax.bar(np.arange(len(max_smd)), max_smd.abs_smd,
           color=[COLORS[t] for t in max_smd.trajectory])
    ax.axhline(0.1, color="#777777", lw=0.75, ls="--")
    ax.set_xticks(np.arange(len(max_smd)))
    ax.set_xticklabels(max_smd.trajectory + " " + max_smd.role,
                       rotation=45, ha="right")
    ax.set_ylabel("Maximum |SMD| after matching")
    clean_ax(ax)
    panel_label(ax, "E", x=-0.23)

    ax = fig.add_subplot(gs[2, 1])
    names = {"any_loop_i": "Any loop", "GGL_i": "GGL", "DGL_i": "DGL"}
    y = np.arange(len(loop))[::-1]
    for yi, (_, row) in zip(y, loop.iterrows()):
        x, lo, hi = row.beta_log2_TPM, row.ci95_low, row.ci95_high
        ax.errorbar(x, yi, xerr=[[x - lo], [hi - x]], fmt="o",
                    color="#555555", capsize=2, lw=0.9, ms=4)
        ax.text(hi + 0.012, yi, f"partial $R^2$={row.partial_R2:.4f}",
                va="center", fontsize=6)
    ax.axvline(0, color="#777777", lw=0.75, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([names[t] for t in loop.term])
    ax.set_xlabel("Adjusted β: log$_2$(TPM + 1)")
    ax.set_xlim(min(-0.17, loop.ci95_low.min() - 0.03),
                max(0.34, loop.ci95_high.max() + 0.12))
    clean_ax(ax)
    panel_label(ax, "F", x=-0.23)

    fig.subplots_adjust(left=0.17, right=0.98, top=0.98, bottom=0.10)
    export(fig, outdir, "Supplementary_Figure_S16_matched_controls_and_loops")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--te", type=Path)
    args = parser.parse_args()
    configure_style()
    fig_s13(args.input, args.output)
    fig_s14(args.input, args.output)
    fig_s15(args.input, args.output, args.te)
    fig_s16(args.input, args.output)


if __name__ == "__main__":
    main()
