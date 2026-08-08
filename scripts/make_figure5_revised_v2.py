#!/usr/bin/env python3
"""Final layout revision for Figure 5."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

import make_figure5 as base
import make_figure5_revised as rev


def draw_phylogeny(ax, thresholds):
    rev.draw_phylogeny(ax, thresholds)
    abbreviated = {
        "Pyrus bretschneideri": "P. bretschneideri",
        "Malus domestica": "M. domestica",
        "Gillenia trifoliata": "G. trifoliata",
        "Prunus persica": "P. persica",
        "Fragaria vesca": "F. vesca",
        "Vitis vinifera": "V. vinifera",
    }
    for artist in list(ax.texts):
        if artist.get_text().startswith("node "):
            artist.remove()
        elif artist.get_text() in abbreviated:
            artist.set_text(abbreviated[artist.get_text()])
            artist.set_fontstyle("italic")
            artist.set_fontsize(6.8)
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
    connector = {
        "arrowstyle": "-",
        "lw": 0.6,
        "color": "#65717B",
        "connectionstyle": "angle3,angleA=0,angleB=90",
    }
    ax.annotate(
        f"node 1: Ks = {mdo:.3f}",
        xy=(0.69, 0.76),
        xytext=(0.42, 0.94),
        arrowprops=connector,
        fontsize=6.7,
        ha="left",
        va="center",
    )
    ax.annotate(
        f"node 2: Ks = {ppe:.3f}",
        xy=(0.39, 0.50),
        xytext=(0.03, 0.59),
        arrowprops=connector,
        fontsize=6.7,
        ha="left",
        va="center",
    )


def draw_synteny_rule(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    species = [
        ("P. bretschneideri", 0.90),
        ("M. domestica", 0.66),
        ("G. trifoliata", 0.42),
        ("P. persica", 0.18),
    ]
    x_positions = [0.35, 0.49, 0.63, 0.77]
    gene_names = ["A", "B", "C", "D"]
    block_colors = ["#A9D0E5", "#E47761", "#C9B8E2", "#AED0AA"]
    for label, ypos in species:
        ax.text(0.02, ypos, label, va="center", ha="left", fontsize=6.8)
        ax.plot(
            [0.29, 0.83], [ypos, ypos], color="#AEB6BD",
            lw=1.15, solid_capstyle="round", zorder=0
        )
        for xpos, gene, color in zip(x_positions, gene_names, block_colors):
            absent = label != "P. bretschneideri" and gene == "B"
            ax.add_patch(
                Rectangle(
                    (xpos - 0.029, ypos - 0.029), 0.058, 0.058,
                    facecolor="white" if absent else color,
                    edgecolor="#E47761" if absent else "white",
                    lw=0.85 if absent else 0.55,
                    linestyle="--" if absent else "-",
                    zorder=2,
                )
            )
            if label == "P. bretschneideri":
                ax.text(
                    xpos, ypos + 0.063, gene,
                    ha="center", va="bottom", fontsize=6.5
                )


def stacked_age_bars(ax, percentages, duplication_type, show_ylabel=False):
    base.stacked_age_bars(
        ax, percentages, duplication_type, show_ylabel=show_ylabel
    )
    ax.set_ylim(0, 114)


def main():
    base.configure()
    thresholds = pd.read_csv(
        base.SRC / "Figure5A_phylogeny_thresholds.tsv", sep="\t"
    )
    percentages = pd.read_csv(
        base.SRC / "Figure5C_trajectory_percentages.tsv", sep="\t"
    )
    ks = pd.read_csv(base.SRC / "Figure5C_Ks.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 8.7))
    outer = fig.add_gridspec(
        2, 2,
        height_ratios=[0.36, 0.64],
        width_ratios=[1.08, 0.92],
        left=0.085, right=0.98, top=0.97, bottom=0.075,
        wspace=0.34, hspace=0.24,
    )
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    draw_phylogeny(ax_a, thresholds)
    draw_synteny_rule(ax_b)

    panel_c = outer[1, :].subgridspec(
        2, 3, height_ratios=[1.05, 0.72], wspace=0.16, hspace=0.34
    )
    upper_axes, lower_axes = [], []
    for column, duplication_type in enumerate(base.TYPES):
        ax_top = fig.add_subplot(panel_c[0, column])
        ax_bottom = fig.add_subplot(panel_c[1, column])
        stacked_age_bars(
            ax_top, percentages, duplication_type, column == 0
        )
        base.density_panel(ax_bottom, ks, duplication_type, column == 0)
        upper_axes.append(ax_top)
        lower_axes.append(ax_bottom)

    base.panel_label(fig, ax_a, "A", x_offset=0.055, y_offset=0.005)
    base.panel_label(fig, ax_b, "B", x_offset=0.020, y_offset=0.005)
    base.panel_label(fig, upper_axes[0], "C", x_offset=0.055, y_offset=0.015)

    trajectory_handles = [
        Patch(facecolor=base.TRAJ_COLORS[item], label=item)
        for item in base.TRAJECTORIES
    ]
    first_box = upper_axes[0].get_position(fig)
    fig.legend(
        handles=trajectory_handles,
        ncol=4,
        loc="upper left",
        bbox_to_anchor=(first_box.x0, first_box.y1 - 0.028),
        bbox_transform=fig.transFigure,
        borderaxespad=0,
        columnspacing=1.0,
        handlelength=1.05,
    )
    age_handles = [
        Line2D(
            [0], [0], color=base.AGE_COLORS[age], lw=1.7,
            label=base.AGE_LABELS[age]
        )
        for age in base.AGES
    ]
    fig.legend(
        handles=age_handles,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.53, lower_axes[0].get_position(fig).y1 + 0.025),
        columnspacing=1.2,
    )
    base.save_bundle(fig, base.OUT / "Figure5_hapA_two_method_consensus")


if __name__ == "__main__":
    main()
