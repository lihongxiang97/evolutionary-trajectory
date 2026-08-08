#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from PIL import Image
from figure3_helpers import (
    CLASSES, SHORT, ROLES, RCOL, CCOL, bh, stars, bracket, draw_box
)

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data" / "main_figures" / "figure3"
OUT = ROOT / "output" / "figure3"
OUT.mkdir(parents=True, exist_ok=True)
TYPES = ["TD", "PD", "TRD"]
TAU_COL = {"P": "parent_tau", "C": "child_tau", "A": "ancestor_tau"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.usetex": False,
})


def draw_tau_grid(fig, grid, data):
    tests, plot_rows = [], []
    axes = {}
    for ri, trajectory in enumerate(CLASSES):
        for ci, dtype in enumerate(TYPES):
            ax = fig.add_subplot(grid[ri, ci])
            axes[ri, ci] = ax
            sub = data[
                (data["two_method_agreement"] == trajectory)
                & (data["duplication_type"] == dtype)
            ]
            values = {
                role: pd.to_numeric(sub[TAU_COL[role]], errors="coerce")
                .dropna()
                .to_numpy()
                for role in ROLES
            }
            draw_box(
                ax,
                [values[role] for role in ROLES],
                [RCOL[role] for role in ROLES],
                ROLES,
                0.54,
            )
            ax.set_ylim(0, 1.34)
            ax.set_yticks(np.arange(0, 1.01, 0.2))
            if ri == 0:
                ax.set_title(dtype, fontsize=9, weight="bold", pad=5)
            if ci == 0:
                # Keep the complete axis title as one editable text object.
                # Unicode tau avoids splitting normal text and mathtext on export.
                ax.set_ylabel("Tissue specificity index (τ)")
                ax.text(
                    -0.38, 0.5, SHORT[trajectory], transform=ax.transAxes,
                    ha="right", va="center", fontsize=9, weight="bold"
                )
            else:
                ax.set_yticklabels([])
            ax.text(
                0.04, 0.95, f"n={len(sub)}", transform=ax.transAxes,
                ha="left", va="top", fontsize=6.2, color="#555555"
            )
            local = []
            for left, right in [("P", "C"), ("P", "A"), ("C", "A")]:
                pair = sub[[TAU_COL[left], TAU_COL[right]]].apply(
                    pd.to_numeric, errors="coerce"
                ).dropna()
                try:
                    statistic, p_value = wilcoxon(pair.iloc[:, 0], pair.iloc[:, 1])
                except ValueError:
                    statistic, p_value = np.nan, 1.0
                row = {
                    "duplication_type": dtype,
                    "trajectory": SHORT[trajectory],
                    "comparison": f"{left} vs {right}",
                    "n": len(pair),
                    "statistic": statistic,
                    "p_raw": p_value,
                }
                tests.append(row)
                local.append(row)
            ax._local_tests = local
            for role in ROLES:
                plot_rows.extend(
                    {
                        "duplication_type": dtype,
                        "trajectory": SHORT[trajectory],
                        "role": role,
                        "tau": value,
                    }
                    for value in values[role]
                )
    q_values = bh([row["p_raw"] for row in tests])
    for row, q_value in zip(tests, q_values):
        row["p_BH"] = q_value
        row["significance"] = stars(q_value)
    for ax in axes.values():
        local = ax._local_tests
        bracket(ax, 0, 1, 1.035, 0.025, local[0]["significance"], 5.2)
        bracket(ax, 1, 2, 1.035, 0.025, local[2]["significance"], 5.2)
        bracket(ax, 0, 2, 1.20, 0.025, local[1]["significance"], 5.2)
    return pd.DataFrame(tests), pd.DataFrame(plot_rows)


def draw_kaks_row(fig, grid, data):
    tests, plot_rows = [], []
    all_kaks = pd.to_numeric(data["Ka_Ks"], errors="coerce").dropna().to_numpy()
    common_upper = 1.6
    for ci, dtype in enumerate(TYPES):
        ax = fig.add_subplot(grid[0, ci])
        subset = data[data["duplication_type"] == dtype]
        values = [
            pd.to_numeric(
                subset.loc[subset["trajectory"] == trajectory, "Ka_Ks"],
                errors="coerce",
            ).dropna().to_numpy()
            for trajectory in CLASSES
        ]
        draw_box(
            ax, values, [CCOL[x] for x in CLASSES],
            [SHORT[x] for x in CLASSES], 0.50
        )
        ax.set_title(dtype, fontsize=9, weight="bold", pad=5)
        if ci == 0:
            ax.set_ylabel("Ka/Ks")
        else:
            ax.set_yticklabels([])
        upper = common_upper
        ax.set_ylim(0, upper)
        for x, vals in enumerate(values):
            ax.text(
                x, 1.5, f"n={len(vals)}",
                ha="center", va="center", fontsize=5.7, color="#555555"
            )
        for trajectory, vals in zip(CLASSES, values):
            plot_rows.extend(
                {
                    "duplication_type": dtype,
                    "trajectory": SHORT[trajectory],
                    "Ka_Ks": value,
                }
                for value in vals
            )
    return pd.DataFrame(tests), pd.DataFrame(plot_rows)


def main():
    tau = pd.read_csv(SRC / "two_method_trajectory_tau.tsv", sep="\t")
    kaks = pd.read_csv(SRC / "two_method_trajectory_kaks.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 10.2))
    outer = fig.add_gridspec(
        2, 1, height_ratios=[4.0, 1.12],
        left=0.17, right=0.985, top=0.965, bottom=0.065, hspace=0.24
    )
    tau_grid = outer[0].subgridspec(4, 3, wspace=0.12, hspace=0.16)
    kaks_grid = outer[1].subgridspec(1, 3, wspace=0.12)
    tau_tests, tau_plot = draw_tau_grid(fig, tau_grid, tau)
    kaks_tests, kaks_plot = draw_kaks_row(fig, kaks_grid, kaks)
    fig.text(0.018, 0.965, "A", fontsize=13, weight="bold", va="top")
    fig.text(
        0.018, outer[1].get_position(fig).y1 + 0.025,
        "B", fontsize=13, weight="bold", va="top"
    )
    base = OUT / "Figure3_hapA_two_method_consensus"
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=300)
    fig.savefig(
        base.with_suffix(".tiff"), dpi=600,
        pil_kwargs={"compression": "tiff_lzw"}
    )
    plt.close(fig)
    tau_tests.to_csv(
        OUT / "Figure3A_tau_pairwise_wilcoxon_BH.tsv", sep="\t", index=False
    )
    tau_plot.to_csv(OUT / "Figure3A_tau_plot_data.tsv", sep="\t", index=False)
    kaks_tests.to_csv(
        OUT / "Figure3B_kaks_pairwise_mannwhitney_BH.tsv",
        sep="\t", index=False
    )
    kaks_plot.to_csv(OUT / "Figure3B_kaks_plot_data.tsv", sep="\t", index=False)
    summary = {
        "two_method_consensus_pairs": len(tau),
        "valid_kaks_pairs": len(kaks),
        "tau_tests": len(tau_tests),
        "kaks_tests": len(kaks_tests),
        "tau_panel_counts": (
            tau.groupby(["two_method_agreement", "duplication_type"])
            .size().to_dict()
        ),
        "kaks_panel_counts": (
            kaks.groupby(["trajectory", "duplication_type"]).size().to_dict()
        ),
    }
    # JSON cannot encode tuple keys.
    summary["tau_panel_counts"] = {
        "|".join(key): value for key, value in summary["tau_panel_counts"].items()
    }
    summary["kaks_panel_counts"] = {
        "|".join(key): value for key, value in summary["kaks_panel_counts"].items()
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





