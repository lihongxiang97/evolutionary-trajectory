#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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
COMP_COLORS = {"A": "#E9C46A", "B": "#4C78A8"}
TAD_COLORS = {"Boundary": "#D65F5F", "Interior": "#66A182"}


def stacked_percent(ax, data, duplication_type, category_column, categories, colors, show_ylabel):
    block = data[data["duplication_type"] == duplication_type]
    bottom = np.zeros(2)
    for category in categories:
        values = []
        for age in AGES:
            hit = block[
                (block["age"] == age) & (block[category_column] == category)
            ]
            values.append(float(hit["percentage"].iloc[0]) if len(hit) else 0)
        ax.bar(
            [0, 1],
            values,
            bottom=bottom,
            width=0.62,
            color=colors[category],
            edgecolor="white",
            linewidth=0.6,
        )
        bottom += np.asarray(values)
    totals = {
        age: int(block.loc[block["age"] == age, "count"].sum()) for age in AGES
    }
    ax.set_xticks(
        [0, 1],
        [f"Young\nn={totals['young']}", f"Old\nn={totals['old']}"],
    )
    ax.set_ylim(0, 100)
    ax.set_ylabel("Genes (%)")


def contact_boxplot(ax, data, tests, duplication_type, show_ylabel):
    block = data[data["duplication_type"] == duplication_type]
    values = [
        block.loc[block["age"] == age, "log2_contact_CPM_plus_1"]
        .dropna()
        .to_numpy()
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
    ax.set_xticks(
        [0, 1],
        [f"Young\nn={len(values[0])}", f"Old\nn={len(values[1])}"],
    )
    ax.set_ylabel("log2(contact CPM + 1)")
    ax.set_ylim(0, 5.35)
    row = tests[
        (tests["duplication_type"] == duplication_type)
        & (tests["metric"] == "TSS-bin contact strength")
    ].iloc[0]
    y = 4.88
    ax.plot([0, 0, 1, 1], [y - 0.08, y, y, y - 0.08], color="#333333", lw=0.6)
    ax.text(
        0.5,
        y + 0.05,
        significance(row["p_BH"]),
        ha="center",
        va="bottom",
        fontsize=7,
    )


def main():
    configure()
    compartment = pd.read_csv(SRC / "Figure7A_compartment.tsv", sep="\t")
    tad = pd.read_csv(SRC / "Figure7B_TAD.tsv", sep="\t")
    contacts = pd.read_csv(SRC / "Figure7C_contacts.tsv", sep="\t")
    tests = pd.read_csv(SRC / "Figure7_tests_BH.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 7.0))
    grid = fig.add_gridspec(
        3,
        3,
        left=0.09,
        right=0.985,
        top=0.935,
        bottom=0.08,
        wspace=0.34,
        hspace=0.48,
    )
    rows = [[], [], []]
    for column, duplication_type in enumerate(TYPES):
        for row in range(3):
            rows[row].append(fig.add_subplot(grid[row, column]))
        rows[0][column].set_title(
            duplication_type, fontsize=9, weight="bold", pad=4
        )
        stacked_percent(
            rows[0][column],
            compartment,
            duplication_type,
            "compartment",
            ["A", "B"],
            COMP_COLORS,
            column == 0,
        )
        stacked_percent(
            rows[1][column],
            tad,
            duplication_type,
            "TAD_location",
            ["Boundary", "Interior"],
            TAD_COLORS,
            column == 0,
        )
        contact_boxplot(
            rows[2][column], contacts, tests, duplication_type, column == 0
        )

    for label, row in zip(["A", "B", "C"], rows):
        panel_label(fig, row[0], label, x_offset=0.055, y_offset=0.012)

    fig.legend(
        handles=[
            Patch(facecolor=COMP_COLORS[category], label=category)
            for category in ["A", "B"]
        ],
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.52, 0.992),
    )
    fig.legend(
        handles=[
            Patch(facecolor=TAD_COLORS[category], label=category)
            for category in ["Boundary", "Interior"]
        ],
        ncol=2,
        loc="center right",
        bbox_to_anchor=(0.985, rows[1][-1].get_position(fig).y1 + 0.025),
    )
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="s",
                color="none",
                markerfacecolor=AGE_COLORS[age],
                markeredgecolor="none",
                markersize=6,
                label=AGE_LABELS[age],
            )
            for age in AGES
        ],
        ncol=2,
        loc="center right",
        bbox_to_anchor=(0.985, rows[2][-1].get_position(fig).y1 + 0.025),
    )
    save_bundle(fig, OUT / "Figure7_hapA")


if __name__ == "__main__":
    main()

