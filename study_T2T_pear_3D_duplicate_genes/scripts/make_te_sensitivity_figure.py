#!/usr/bin/env python3
"""Create a publication-grade TE-overlap and TE-adjusted TRD sensitivity figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


YOUNG = "#D95F02"
OLD = "#4C78A8"
BASELINE = "#777777"
TE_ADJUSTED = "#7B3294"


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def clean(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.7, length=3)


def errorbar(ax, row, y, color, label=None) -> None:
    x = row.effect_young_vs_old
    ax.errorbar(
        x,
        y,
        xerr=[[x - row.ci95_low], [row.ci95_high - x]],
        fmt="o",
        ms=4,
        lw=0.9,
        capsize=2,
        color=color,
        label=label,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    style()

    te = pd.read_csv(args.results / "TE_overlap_summary.tsv", sep="\t")
    te_models = pd.read_csv(
        args.results / "TE_adjusted_TRD_age_models.tsv", sep="\t"
    )
    baseline = pd.read_csv(
        args.results / "adjusted_young_old_chromatin_models.tsv", sep="\t"
    )
    te = te[
        te.duplication_type.eq("TRD") & te.feature_type.eq("binary_overlap")
    ].set_index("feature")
    te_models = te_models[
        te_models.sensitivity.eq("complete_ORF_expressed")
    ].set_index("outcome")
    baseline = baseline[
        baseline.duplication_type.eq("TRD")
        & baseline.sensitivity.eq("complete_ORF_expressed")
    ].set_index("outcome")

    fig, axes = plt.subplots(
        1, 3, figsize=(7.2, 2.55), gridspec_kw={"width_ratios": [1.05, 1.35, 0.85]}
    )

    ax = axes[0]
    features = ["gene_body_TE", "promoter_TE"]
    labels = ["Gene body", "2-kb promoter"]
    x = np.arange(2)
    width = 0.34
    young = [100 * te.loc[f, "young_fraction"] for f in features]
    old = [100 * te.loc[f, "old_fraction"] for f in features]
    ax.bar(x - width / 2, young, width, color=YOUNG, label="Young")
    ax.bar(x + width / 2, old, width, color=OLD, label="Old")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("TRD genes overlapping a TE (%)")
    ymax = max(young + old) * 1.24
    ax.set_ylim(0, ymax)
    for i, f in enumerate(features):
        row = te.loc[f]
        ax.text(
            i,
            ymax * 0.94,
            f"OR={row.effect:.2f}\nBH P={row.p_BH:.2g}",
            ha="center",
            va="top",
            fontsize=6.2,
        )
    clean(ax)
    ax.text(-0.22, 1.05, "A", transform=ax.transAxes, fontweight="bold", fontsize=10)

    ax = axes[1]
    outcomes = [
        ("log2_TPM_plus_1", "Expression"),
        ("mean_ATAC_TSS_2kb", "ATAC"),
        ("PC1", "PC1"),
        ("log2_contact_CPM_plus_1", "Contact"),
    ]
    y = np.arange(len(outcomes))[::-1]
    for yi, (outcome, _) in zip(y, outcomes):
        errorbar(ax, baseline.loc[outcome], yi + 0.12, BASELINE)
        errorbar(ax, te_models.loc[outcome], yi - 0.12, TE_ADJUSTED)
    ax.axvline(0, color="#999999", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([label for _, label in outcomes])
    ax.set_xlabel("Adjusted Young–Old coefficient (95% CI)")
    clean(ax)
    ax.text(-0.25, 1.05, "B", transform=ax.transAxes, fontweight="bold", fontsize=10)

    ax = axes[2]
    outcomes = [("compartment_A", "A compartment"), ("TAD_boundary", "TAD boundary")]
    y = np.arange(len(outcomes))[::-1]
    for yi, (outcome, _) in zip(y, outcomes):
        errorbar(ax, baseline.loc[outcome], yi + 0.12, BASELINE)
        errorbar(ax, te_models.loc[outcome], yi - 0.12, TE_ADJUSTED)
    ax.axvline(1, color="#999999", lw=0.7, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([label for _, label in outcomes])
    ax.set_xlabel("Adjusted odds ratio (95% CI)")
    clean(ax)
    ax.text(-0.30, 1.05, "C", transform=ax.transAxes, fontweight="bold", fontsize=10)

    handles = [
        Line2D([0], [0], marker="o", color=BASELINE, lw=0, label="Genomic covariates"),
        Line2D([0], [0], marker="o", color=TE_ADJUSTED, lw=0, label="+ TE coverage"),
    ]
    fig.legend(
        handles=handles,
        frameon=False,
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.66, -0.01),
    )
    age_handles = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor={"Young": YOUNG, "Old": OLD}[age],
               markeredgecolor="none", markersize=7, label=age)
        for age in ["Young", "Old"]
    ]
    fig.legend(
        handles=age_handles,
        frameon=False,
        ncol=2,
        loc="upper left",
        bbox_to_anchor=(0.08, 0.985),
    )
    fig.subplots_adjust(left=0.085, right=0.985, top=0.82, bottom=0.24, wspace=0.62)

    args.output.mkdir(parents=True, exist_ok=True)
    stem = args.output / "Supplementary_Figure_S17_TE_sensitivity"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=300,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
