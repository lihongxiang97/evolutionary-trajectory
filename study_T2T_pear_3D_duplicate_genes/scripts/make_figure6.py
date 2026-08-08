#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figure_utils import (
    AGE_COLORS,
    AGE_LABELS,
    configure,
    panel_label,
    save_bundle,
    significance,
)


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "source"
OUT = ROOT / "output"
TYPES = ["TD", "PD", "TRD"]
AGES = ["young", "old"]


def expression_boxplot(ax, data, tests, duplication_type, show_ylabel=False):
    block = data[data["duplication_type"] == duplication_type]
    values = [
        block.loc[block["age"] == age, "log2_TPM_plus_1"].dropna().to_numpy()
        for age in AGES
    ]
    bp = ax.boxplot(
        values,
        positions=[0, 1],
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "white", "lw": 0.9},
        whiskerprops={"lw": 0.65},
        capprops={"lw": 0.65},
        boxprops={"lw": 0.65},
    )
    for patch, age in zip(bp["boxes"], AGES):
        patch.set_facecolor(AGE_COLORS[age])
        patch.set_alpha(0.92)
    ax.set_xticks(
        [0, 1],
        [f"Young\nn={len(values[0])}", f"Old\nn={len(values[1])}"],
    )
    ax.set_title(duplication_type, fontsize=9, weight="bold", pad=4)
    ax.set_ylabel("log2(TPM + 1)")

    row = tests[tests["duplication_type"] == duplication_type].iloc[0]
    ymax = max(np.quantile(values[0], 0.98), np.quantile(values[1], 0.98))
    ax.set_ylim(0, ymax * 1.22)
    y = ymax * 1.08
    ax.plot([0, 0, 1, 1], [y * 0.98, y, y, y * 0.98], color="#333333", lw=0.6)
    ax.text(
        0.5,
        y * 1.01,
        significance(row["p_BH"]),
        ha="center",
        va="bottom",
        fontsize=7,
    )


def atac_profile(ax, profile, duplication_type, show_ylabel=False):
    block = profile[profile["duplication_type"] == duplication_type]
    for age in AGES:
        sub = block[block["age"] == age].sort_values("profile_bin")
        x = sub["profile_bin"].to_numpy()
        ax.plot(
            x,
            sub["mean"],
            color=AGE_COLORS[age],
            lw=1.35,
            label=AGE_LABELS[age],
        )
        ax.fill_between(
            x,
            sub["ci95_low"],
            sub["ci95_high"],
            color=AGE_COLORS[age],
            alpha=0.15,
            lw=0,
        )
    ax.axvline(19.5, color="#777777", ls="--", lw=0.55)
    ax.axvline(59.5, color="#777777", ls="--", lw=0.55)
    ax.set_xlim(0, 79)
    ax.set_xticks([0, 19.5, 59.5, 79], ["−2 kb", "TSS", "TES", "+2 kb"])
    ax.set_xlabel("Scaled gene region")
    ax.set_ylabel("ATAC signal (RPGC)")
    counts = {
        age: int(
            block.loc[block["age"] == age, "count"].dropna().max()
        )
        for age in AGES
    }
    ax.text(
        0.98,
        1.03,
        f"n={counts['young']}/{counts['old']}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.8,
        color="#555555",
    )


def main():
    configure()
    expression = pd.read_csv(SRC / "Figure6A_expression.tsv", sep="\t")
    expression_tests = pd.read_csv(
        SRC / "Figure6A_expression_tests_BH.tsv", sep="\t"
    )
    profile = pd.read_csv(SRC / "Figure6B_ATAC_profiles.tsv.gz", sep="\t")

    fig = plt.figure(figsize=(7.2, 4.65))
    grid = fig.add_gridspec(
        2,
        3,
        left=0.09,
        right=0.985,
        top=0.91,
        bottom=0.12,
        hspace=0.58,
        wspace=0.34,
    )
    top_axes = []
    bottom_axes = []
    for column, duplication_type in enumerate(TYPES):
        ax_a = fig.add_subplot(grid[0, column])
        ax_b = fig.add_subplot(grid[1, column])
        expression_boxplot(
            ax_a, expression, expression_tests, duplication_type, column == 0
        )
        atac_profile(ax_b, profile, duplication_type, column == 0)
        top_axes.append(ax_a)
        bottom_axes.append(ax_b)

    panel_label(fig, top_axes[0], "A", x_offset=0.055, y_offset=0.012)
    panel_label(fig, bottom_axes[0], "B", x_offset=0.055, y_offset=0.012)
    handles = [
        Line2D(
            [0],
            [0],
            color=AGE_COLORS[age],
            lw=1.7,
            label=AGE_LABELS[age],
        )
        for age in AGES
    ]
    fig.legend(
        handles=handles,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.995),
        columnspacing=1.2,
    )
    save_bundle(fig, OUT / "Figure6_hapA")


if __name__ == "__main__":
    main()

