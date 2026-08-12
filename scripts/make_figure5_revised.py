#!/usr/bin/env python3
"""Revised Figure 5 with a polished phylogeny and uncluttered synteny panel."""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import make_figure5 as base


SPECIES_COLORS = {
    "Pyrus bretschneideri": "#D65F5F",
    "Malus domestica": "#E89B3C",
    "Gillenia trifoliata": "#7A9E5B",
    "Prunus persica": "#4C9F9A",
    "Fragaria vesca": "#5C78B5",
    "Vitis vinifera": "#8C6BB1",
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
    y = dict(species)
    x_tip = 0.91
    branch = "#64707A"

    node_x = 0.69
    clade_y = (y["Pyrus bretschneideri"] + y["Malus domestica"]) / 2
    for name in ("Pyrus bretschneideri", "Malus domestica"):
        ax.plot([node_x, x_tip], [y[name], y[name]], color=branch, lw=1.15)
    ax.plot(
        [node_x, node_x],
        [y["Malus domestica"], y["Pyrus bretschneideri"]],
        color=branch,
        lw=1.15,
    )
    ax.scatter(
        node_x, clade_y, s=17, facecolor="white",
        edgecolor=branch, lw=0.8, zorder=4
    )
    clade_x = node_x

    for name, next_x in [
        ("Gillenia trifoliata", 0.52),
        ("Prunus persica", 0.39),
        ("Fragaria vesca", 0.26),
        ("Vitis vinifera", 0.13),
    ]:
        tip_y = y[name]
        old_clade_y = clade_y
        ax.plot([next_x, x_tip], [tip_y, tip_y], color=branch, lw=1.15)
        ax.plot([next_x, clade_x], [old_clade_y, old_clade_y], color=branch, lw=1.15)
        ax.plot([next_x, next_x], [tip_y, old_clade_y], color=branch, lw=1.15)
        clade_y = (tip_y + old_clade_y) / 2
        clade_x = next_x
        ax.scatter(
            clade_x, clade_y, s=14, facecolor="white",
            edgecolor=branch, lw=0.75, zorder=4
        )

    for name, ypos in species:
        color = SPECIES_COLORS[name]
        ax.scatter(
            x_tip, ypos, s=29, facecolor=color,
            edgecolor="white", lw=0.6, zorder=5
        )
        ax.text(
            x_tip + 0.026, ypos, name, va="center", ha="left",
            fontsize=7.1, color=color, weight="medium"
        )

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
        xy=(0.69, 0.76),
        xytext=(0.29, 0.87),
        arrowprops={"arrowstyle": "-", "lw": 0.65, "color": "#65717B"},
        fontsize=6.7,
        ha="left",
    )
    ax.annotate(
        f"node 2: Ks = {ppe:.3f}",
        xy=(0.39, 0.50),
        xytext=(0.05, 0.47),
        arrowprops={"arrowstyle": "-", "lw": 0.65, "color": "#65717B"},
        fontsize=6.7,
        ha="left",
    )
    ax.scatter(
        0.615, 0.76, marker="*", s=88, facecolor="#F2C14E",
        edgecolor="#9B6A00", lw=0.65, zorder=6
    )
    ax.text(
        0.585, 0.76, "recent WGD", fontsize=6.4, color="#8A6100",
        rotation=90, ha="right", va="center"
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
                    xpos, ypos + 0.066, gene,
                    ha="center", va="bottom", fontsize=6.5
                )


if __name__ == "__main__":
    base.draw_phylogeny = draw_phylogeny
    base.draw_synteny_rule = draw_synteny_rule
    base.main()
