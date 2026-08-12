#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, Patch, Rectangle
from scipy.stats import gaussian_kde

from figure_utils import AGE_COLORS, AGE_LABELS, configure, panel_label, save_bundle


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data" / "main_figures" / "figures5_7"
OUT = ROOT / "output" / "figure5"
OUT.mkdir(parents=True, exist_ok=True)

TYPES = ["TD", "PD", "TRD"]
AGES = ["young", "old"]
TRAJECTORIES = ["Con", "Neo(C)", "Neo(P)", "Spe"]
TRAJ_COLORS = {
    "Con": "#4C78A8",
    "Neo(C)": "#F2A541",
    "Neo(P)": "#9C6ADE",
    "Spe": "#D65F5F",
}


def draw_phylogeny(ax, thresholds):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    species = [
        ("Pyrus bretschneideri", 0.84),
        ("Malus domestica", 0.68),
        ("Gillenia trifoliata", 0.52),
        ("Prunus persica", 0.36),
        ("Fragaria vesca", 0.20),
        ("Vitis vinifera", 0.06),
    ]
    x_tip = 0.93
    y = {name: value for name, value in species}

    # Pyrus and Malus are sisters; successive nodes add the remaining taxa.
    node_x = 0.69
    clade_y = (y["Pyrus bretschneideri"] + y["Malus domestica"]) / 2
    ax.plot([node_x, x_tip], [y["Pyrus bretschneideri"]] * 2, color="#4D4D4D", lw=1.1)
    ax.plot([node_x, x_tip], [y["Malus domestica"]] * 2, color="#4D4D4D", lw=1.1)
    ax.plot([node_x, node_x], [y["Malus domestica"], y["Pyrus bretschneideri"]], color="#4D4D4D", lw=1.1)
    clade_x = node_x

    for name, next_x in [
        ("Gillenia trifoliata", 0.52),
        ("Prunus persica", 0.39),
        ("Fragaria vesca", 0.26),
        ("Vitis vinifera", 0.13),
    ]:
        tip_y = y[name]
        ax.plot([next_x, x_tip], [tip_y, tip_y], color="#4D4D4D", lw=1.1)
        ax.plot([next_x, clade_x], [clade_y, clade_y], color="#4D4D4D", lw=1.1)
        ax.plot([next_x, next_x], [tip_y, clade_y], color="#4D4D4D", lw=1.1)
        clade_y = (tip_y + clade_y) / 2
        clade_x = next_x

    for name, ypos in species:
        ax.text(x_tip + 0.02, ypos, name, va="center", ha="left", fontsize=7.2)

    mdo = float(
        thresholds.loc[
            thresholds["comparison"].eq("Pbr_hapA-Mdo"),
            "mean_Ks_primary_threshold",
        ].iloc[0]
    )
    ppe = float(
        thresholds.loc[
            thresholds["comparison"].eq("Pbr_hapA-Ppe"),
            "mean_Ks_primary_threshold",
        ].iloc[0]
    )
    ax.annotate(
        f"node 1: Ks = {mdo:.3f}",
        xy=(0.69, (0.84 + 0.68) / 2),
        xytext=(0.29, 0.83),
        arrowprops={"arrowstyle": "-", "lw": 0.7, "color": "#555555"},
        fontsize=6.8,
        ha="left",
    )
    ax.annotate(
        f"node 2: Ks = {ppe:.3f}",
        xy=(0.39, 0.50),
        xytext=(0.05, 0.48),
        arrowprops={"arrowstyle": "-", "lw": 0.7, "color": "#555555"},
        fontsize=6.8,
        ha="left",
    )
    ax.text(
        0.60,
        0.75,
        "recent WGD",
        fontsize=6.5,
        color="#777777",
        rotation=90,
        ha="center",
        va="center",
    )


def draw_synteny_rule(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    species = [
        ("P. bretschneideri", 0.82),
        ("M. domestica", 0.62),
        ("G. trifoliata", 0.42),
        ("P. persica", 0.22),
    ]
    x_positions = [0.34, 0.48, 0.62, 0.76]
    gene_names = ["A", "B", "C", "D"]
    block_colors = ["#B9D6E8", "#E07A5F", "#D4C4E8", "#BFD8BD"]

    for label, ypos in species:
        ax.text(0.03, ypos, label, va="center", ha="left", fontsize=7)
        ax.plot([0.29, 0.83], [ypos, ypos], color="#B8B8B8", lw=5, solid_capstyle="round")
        for idx, (xpos, gene, color) in enumerate(
            zip(x_positions, gene_names, block_colors)
        ):
            absent = label != "P. bretschneideri" and gene == "B"
            if absent:
                ax.add_patch(
                    Rectangle(
                        (xpos - 0.027, ypos - 0.027),
                        0.054,
                        0.054,
                        facecolor="white",
                        edgecolor="#E07A5F",
                        lw=0.8,
                        linestyle="--",
                    )
                )
            else:
                ax.add_patch(
                    Rectangle(
                        (xpos - 0.027, ypos - 0.027),
                        0.054,
                        0.054,
                        facecolor=color,
                        edgecolor="white",
                        lw=0.5,
                    )
                )
            if label == "P. bretschneideri":
                ax.text(xpos, ypos + 0.065, gene, ha="center", va="bottom", fontsize=6.5)

    ax.add_patch(
        FancyBboxPatch(
            (0.29, 0.035),
            0.54,
            0.09,
            boxstyle="round,pad=0.012,rounding_size=0.015",
            facecolor="#F6F6F6",
            edgecolor="#B5B5B5",
            lw=0.7,
        )
    )
    ax.text(
        0.56,
        0.08,
        "Young gene: absent from all three syntenic outgroups\nand pair Ks < node 1",
        ha="center",
        va="center",
        fontsize=6.3,
    )
    


def stacked_age_bars(ax, percentages, duplication_type, show_ylabel=False):
    x = np.array([0, 1])
    bottom = np.zeros(2)
    for trajectory in TRAJECTORIES:
        values = []
        for age in AGES:
            hit = percentages[
                (percentages["duplication_type"] == duplication_type)
                & (percentages["pair_age"] == age)
                & (percentages["trajectory"] == trajectory)
            ]
            values.append(float(hit["percentage"].iloc[0]) if len(hit) else 0)
        ax.bar(
            x,
            values,
            bottom=bottom,
            width=0.66,
            color=TRAJ_COLORS[trajectory],
            edgecolor="white",
            linewidth=0.6,
        )
        bottom += np.asarray(values)
    ax.set_ylim(0, 100)
    ax.set_xticks(x, ["Young", "Old"])
    ax.set_title(duplication_type, fontsize=9, weight="bold", pad=4)
    if show_ylabel:
        ax.set_ylabel("Pairs (%)")
    else:
        ax.set_yticklabels([])


def density_panel(ax, ks, duplication_type, show_ylabel=False):
    block = ks[ks["duplication_type"] == duplication_type]
    upper = max(0.15, float(block["Ks"].quantile(0.995)) * 1.05)
    x = np.linspace(0, upper, 300)
    for age in AGES:
        values = block.loc[block["pair_age"] == age, "Ks"].dropna().to_numpy()
        if len(values) > 2 and np.std(values) > 0:
            kde = gaussian_kde(values)
            kde.set_bandwidth(kde.factor * 0.9)
            y = kde(x)
            ax.plot(x, y, color=AGE_COLORS[age], lw=1.45)
            ax.fill_between(x, 0, y, color=AGE_COLORS[age], alpha=0.13, lw=0)
        ax.text(
            0.98,
            0.88 if age == "young" else 0.75,
            f"{AGE_LABELS[age]} n={len(values)}",
            color=AGE_COLORS[age],
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=5.8,
        )
    ax.set_xlim(0, upper)
    ax.set_xlabel("Ks")
    if show_ylabel:
        ax.set_ylabel("Density")
    else:
        ax.set_yticklabels([])


def main():
    configure()
    thresholds = pd.read_csv(SRC / "Figure5A_phylogeny_thresholds.tsv", sep="\t")
    percentages = pd.read_csv(
        SRC / "Figure5C_trajectory_percentages.tsv", sep="\t"
    )
    ks = pd.read_csv(SRC / "Figure5C_Ks.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 8.7))
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.36, 0.64],
        width_ratios=[1.18, 0.82],
        left=0.085,
        right=0.98,
        top=0.97,
        bottom=0.075,
        wspace=0.28,
        hspace=0.24,
    )
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    draw_phylogeny(ax_a, thresholds)
    draw_synteny_rule(ax_b)

    panel_c = outer[1, :].subgridspec(
        2, 3, height_ratios=[1.05, 0.72], wspace=0.16, hspace=0.34
    )
    upper_axes = []
    lower_axes = []
    for column, duplication_type in enumerate(TYPES):
        ax_top = fig.add_subplot(panel_c[0, column])
        ax_bottom = fig.add_subplot(panel_c[1, column])
        stacked_age_bars(ax_top, percentages, duplication_type, column == 0)
        density_panel(ax_bottom, ks, duplication_type, column == 0)
        upper_axes.append(ax_top)
        lower_axes.append(ax_bottom)

    panel_label(fig, ax_a, "A", x_offset=0.055, y_offset=0.005)
    panel_label(fig, ax_b, "B", x_offset=0.045, y_offset=0.005)
    panel_label(fig, upper_axes[0], "C", x_offset=0.055, y_offset=0.015)

    trajectory_handles = [
        Patch(facecolor=TRAJ_COLORS[item], label=item) for item in TRAJECTORIES
    ]
    fig.legend(
        handles=trajectory_handles,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.53, upper_axes[0].get_position(fig).y1 + 0.035),
        columnspacing=1.1,
        handlelength=1.1,
    )
    age_handles = [
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
        handles=age_handles,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.53, lower_axes[0].get_position(fig).y1 + 0.025),
        columnspacing=1.2,
    )

    save_bundle(fig, OUT / "Figure5_hapA_two_method_consensus")


if __name__ == "__main__":
    main()
