#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import wilcoxon, mannwhitneyu
from PIL import Image

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "source"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)

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
ROLES = ["P", "C", "A"]
RCOL = {"P": "#4472C4", "C": "#E07A5F", "A": "#595959"}
CCOL = {
    "Conservation": "#4C78A8",
    "Neofunctionalization(Child)": "#F2A541",
    "Neofunctionalization(Parent)": "#D98C32",
    "Specialization": "#D65F5F",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate(
        (ranked * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1]
    result = np.empty_like(q)
    result[order] = np.minimum(q, 1)
    return result


def stars(p):
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def bracket(ax, x1, x2, y, h, label, size=6):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="#333333", lw=0.6)
    ax.text((x1 + x2) / 2, y + h, label, ha="center", va="bottom", fontsize=size)


def draw_box(ax, values, colors, labels, widths=0.56):
    box = ax.boxplot(
        values,
        positions=np.arange(len(values)),
        widths=widths,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "white", "lw": 1.0},
        whiskerprops={"lw": 0.65},
        capprops={"lw": 0.65},
        boxprops={"lw": 0.65},
    )
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.92)
    ax.set_xticks(np.arange(len(labels)), labels)
    return box


def panel_a(fig, grid, data):
    tests = []
    plot_rows = []
    axes = []
    for index, trajectory in enumerate(CLASSES):
        ax = fig.add_subplot(grid[index, 0])
        axes.append(ax)
        sub = data[data["two_method_agreement"] == trajectory]
        values = {
            "P": pd.to_numeric(sub["parent_tau"], errors="coerce").dropna().to_numpy(),
            "C": pd.to_numeric(sub["child_tau"], errors="coerce").dropna().to_numpy(),
            "A": pd.to_numeric(sub["ancestor_tau"], errors="coerce").dropna().to_numpy(),
        }
        draw_box(ax, [values[r] for r in ROLES], [RCOL[r] for r in ROLES], ROLES)
        ax.text(-0.12, 0.5, SHORT[trajectory], transform=ax.transAxes, ha="right", va="center", fontsize=9, weight="bold")
        ax.set_ylim(0, 1.33)
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.set_ylabel(r"Tissue specificity index ($\tau$)")
        ax.text(
            0.04,
            0.96,
            f"n={len(sub)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=6.5,
            color="#555555",
        )
        comparisons = [("P", "C"), ("P", "A"), ("C", "A")]
        local = []
        for left, right in comparisons:
            pair = sub[[{"P": "parent_tau", "C": "child_tau", "A": "ancestor_tau"}[left],
                        {"P": "parent_tau", "C": "child_tau", "A": "ancestor_tau"}[right]]].apply(
                pd.to_numeric, errors="coerce"
            ).dropna()
            try:
                statistic, p_value = wilcoxon(pair.iloc[:, 0], pair.iloc[:, 1])
            except ValueError:
                statistic, p_value = np.nan, 1.0
            row = {
                "trajectory": SHORT[trajectory],
                "comparison": f"{left} vs {right}",
                "n": len(pair),
                "statistic": statistic,
                "p_raw": p_value,
            }
            tests.append(row)
            local.append(row)
        for role in ROLES:
            plot_rows.extend(
                {
                    "trajectory": SHORT[trajectory],
                    "role": role,
                    "tau": value,
                }
                for value in values[role]
            )
        ax._local_tests = local
    q_values = bh([row["p_raw"] for row in tests])
    for row, q_value in zip(tests, q_values):
        row["p_BH"] = q_value
        row["significance"] = stars(q_value)
    for ax in axes:
        local = ax._local_tests
        bracket(ax, 0, 1, 1.04, 0.025, local[0]["significance"])
        bracket(ax, 1, 2, 1.04, 0.025, local[2]["significance"])
        bracket(ax, 0, 2, 1.20, 0.025, local[1]["significance"])
    return pd.DataFrame(tests), pd.DataFrame(plot_rows)


def panel_b(ax, data):
    values = [
        pd.to_numeric(
            data.loc[data["trajectory"] == trajectory, "Ka_Ks"], errors="coerce"
        ).dropna().to_numpy()
        for trajectory in CLASSES
    ]
    draw_box(ax, values, [CCOL[x] for x in CLASSES], [SHORT[x] for x in CLASSES], 0.52)
    ax.set_ylabel("Ka/Ks")
    pooled = np.concatenate(values)
    upper = min(max(np.quantile(pooled, 0.99) * 1.35, 0.6), 2.0)
    ax.set_ylim(0, upper)
    tests = []
    for left in range(4):
        for right in range(left + 1, 4):
            statistic, p_value = mannwhitneyu(
                values[left], values[right], alternative="two-sided"
            )
            tests.append(
                {
                    "comparison": f"{SHORT[CLASSES[left]]} vs {SHORT[CLASSES[right]]}",
                    "n_left": len(values[left]),
                    "n_right": len(values[right]),
                    "statistic": statistic,
                    "p_raw": p_value,
                }
            )
    q_values = bh([row["p_raw"] for row in tests])
    for row, q_value in zip(tests, q_values):
        row["p_BH"] = q_value
        row["significance"] = stars(q_value)
    # Six comparisons are reported in the source table; the plot shows the
    # biologically ordered adjacent contrasts to avoid an unreadable bracket stack.
    lookup = {row["comparison"]: row for row in tests}
    y0 = upper * 0.77
    for level, (x1, x2, key) in enumerate(
        [
            (0, 1, "Con vs Neo(C)"),
            (1, 2, "Neo(C) vs Neo(P)"),
            (2, 3, "Neo(P) vs Spe"),
        ]
    ):
        bracket(
            ax,
            x1,
            x2,
            y0 + level * upper * 0.075,
            upper * 0.018,
            lookup[key]["significance"],
            6.5,
        )
    for x, vals in enumerate(values):
        ax.text(
            x,
            upper * 0.025,
            f"n={len(vals)}",
            ha="center",
            va="bottom",
            fontsize=6.5,
            color="#555555",
        )
    plot_rows = []
    for trajectory, vals in zip(CLASSES, values):
        plot_rows.extend(
            {"trajectory": SHORT[trajectory], "Ka_Ks": value} for value in vals
        )
    return pd.DataFrame(tests), pd.DataFrame(plot_rows)


def main():
    tau = pd.read_csv(SRC / "two_method_trajectory_tau.tsv", sep="\t")
    kaks = pd.read_csv(SRC / "two_method_trajectory_kaks.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 10.9))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[4.0, 1.35],
        left=0.18,
        right=0.985,
        top=0.965,
        bottom=0.065,
        hspace=0.22,
    )
    panel_a_grid = outer[0].subgridspec(4, 1, hspace=0.16)
    tau_tests, tau_plot = panel_a(fig, panel_a_grid, tau)
    ax_b = fig.add_subplot(outer[1])
    kaks_tests, kaks_plot = panel_b(ax_b, kaks)
    fig.text(0.018, 0.965, "A", fontsize=13, weight="bold", va="top")
    fig.text(
        0.018,
        outer[1].get_position(fig).y1 + 0.025,
        "B",
        fontsize=13,
        weight="bold",
        va="top",
    )
    base = OUT / "Figure3_hapA_two_method_consensus"
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=300)
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    tau_tests.to_csv(
        OUT / "Figure3A_tau_pairwise_wilcoxon_BH.tsv", sep="\t", index=False
    )
    tau_plot.to_csv(OUT / "Figure3A_tau_plot_data.tsv", sep="\t", index=False)
    kaks_tests.to_csv(
        OUT / "Figure3B_kaks_pairwise_mannwhitney_BH.tsv", sep="\t", index=False
    )
    kaks_plot.to_csv(OUT / "Figure3B_kaks_plot_data.tsv", sep="\t", index=False)
    summary = {
        "two_method_consensus_pairs": len(tau),
        "valid_kaks_pairs": len(kaks),
        "trajectory_counts": tau["two_method_agreement"].value_counts().to_dict(),
        "kaks_counts": kaks["trajectory"].value_counts().to_dict(),
        "tau_tests": len(tau_tests),
        "kaks_tests": len(kaks_tests),
    }
    (OUT / "Figure3_analysis_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    image = Image.open(base.with_suffix(".tiff"))
    assert image.info["dpi"][0] >= 599
    assert "<text" in base.with_suffix(".svg").read_text(encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()


