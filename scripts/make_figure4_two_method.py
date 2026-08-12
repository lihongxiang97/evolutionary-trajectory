#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import (
    wilcoxon,
    mannwhitneyu,
    fisher_exact,
    binomtest,
    chi2_contingency,
    chi2,
)
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data" / "main_figures" / "figure4"
OUT = ROOT / "output" / "figure4"
OUT.mkdir(parents=True, exist_ok=True)

TRAJECTORIES = ["Con", "Neo(C)", "Neo(P)", "Spe"]
ROLES = ["P", "C", "S"]
ROLE_LABEL = {"P": "Parent", "C": "Child", "S": "Single"}
RCOL = {"P": "#4472C4", "C": "#E07A5F", "S": "#595959"}
COMP_COL = {"A": "#E9C46A", "B": "#4C78A8"}
TAD_COL = {"Boundary": "#D65F5F", "Interior": "#66A182", "Outside": "#D9D9D9"}
LOOP_COL = {"Any loop": "#4C78A8", "GGL": "#66A182", "DGL": "#D65F5F"}
ATAC_YLIM = {
    "Con": (0.8, 1.6),
    "Neo(C)": (0.8, 1.7),
    "Neo(P)": (0.7, 1.6),
    "Spe": (0.8, 1.55),
}
ATAC_YTICKS = {
    "Con": [0.8, 1.0, 1.2, 1.4, 1.6],
    "Neo(C)": [0.8, 1.0, 1.2, 1.4, 1.6],
    "Neo(P)": [0.8, 1.0, 1.2, 1.4, 1.6],
    "Spe": [0.8, 1.0, 1.2, 1.4],
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.2,
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


def finish_tests(rows):
    if not rows:
        return pd.DataFrame()
    q_values = bh([row["p_raw"] for row in rows])
    for row, q_value in zip(rows, q_values):
        row["p_BH"] = q_value
    return pd.DataFrame(rows)


def paired_binary_test(left, right):
    table = pd.crosstab(left.astype(bool), right.astype(bool)).reindex(
        index=[False, True], columns=[False, True], fill_value=0
    )
    b = int(table.loc[False, True])
    c = int(table.loc[True, False])
    discordant = b + c
    p_value = binomtest(min(b, c), discordant, 0.5).pvalue if discordant else 1.0
    return discordant, p_value


def paired_multicategory_bowker(left, right, categories):
    table = pd.crosstab(left, right).reindex(
        index=categories, columns=categories, fill_value=0
    )
    statistic = 0.0
    df = 0
    for i in range(len(categories)):
        for j in range(i + 1, len(categories)):
            nij = table.iloc[i, j]
            nji = table.iloc[j, i]
            if nij + nji:
                statistic += (nij - nji) ** 2 / (nij + nji)
                df += 1
    p_value = chi2.sf(statistic, df) if df else 1.0
    return statistic, df, p_value


def plot_atac(ax, profile, trajectory, show_ylabel=False):
    block = profile[profile["trajectory"] == trajectory]
    for role in ROLES:
        sub = block[block["role"] == role].sort_values("position")
        x = sub["position"].to_numpy() / 1000
        mean = sub["mean"].to_numpy()
        low = sub["ci95_low"].to_numpy()
        high = sub["ci95_high"].to_numpy()
        ax.plot(x, mean, color=RCOL[role], lw=1.25)
        ax.fill_between(x, low, high, color=RCOL[role], alpha=0.15, lw=0)
    ax.axvline(0, color="#777777", lw=0.55, ls="--")
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(*ATAC_YLIM[trajectory])
    ax.set_yticks(ATAC_YTICKS[trajectory])
    ax.set_xticks([-2.5, 0, 2.5], ["-2.5 kb", "TSS", "+2.5 kb"])
    ax.set_xlabel("Distance from TSS")
    ax.set_ylabel("ATAC signal (RPGC)", labelpad=2)


def plot_compartment(ax, block, show_ylabel=False):
    bottom = np.zeros(3)
    for category in ["A", "B"]:
        values = []
        for role in ROLES:
            role_data = block[block["role"] == role]["compartment"]
            valid = role_data[role_data.isin(["A", "B"])]
            values.append(100 * (valid == category).mean())
        ax.bar(
            [ROLE_LABEL[role] for role in ROLES], values, bottom=bottom, width=0.62,
            color=COMP_COL[category], edgecolor="white", linewidth=0.6
        )
        bottom += np.asarray(values)
    ax.set_ylim(0, 100)
    if show_ylabel:
        ax.set_ylabel("Genes (%)")
    else:
        ax.set_yticklabels([])


def plot_tad(ax, block, show_ylabel=False):
    bottom = np.zeros(3)
    for category in ["Boundary", "Interior", "Outside"]:
        values = []
        for role in ROLES:
            role_data = block[block["role"] == role]["TAD_location"]
            values.append(100 * (role_data == category).mean())
        ax.bar(
            [ROLE_LABEL[role] for role in ROLES], values, bottom=bottom, width=0.62,
            color=TAD_COL[category], edgecolor="white", linewidth=0.6
        )
        bottom += np.asarray(values)
    ax.set_ylim(0, 100)
    if show_ylabel:
        ax.set_ylabel("Genes (%)")
    else:
        ax.set_yticklabels([])


def plot_loops(ax, block, show_ylabel=False):
    x = np.arange(3)
    width = 0.22
    for offset, (label, column) in enumerate(
        [("Any loop", "any_loop"), ("GGL", "GGL"), ("DGL", "DGL")]
    ):
        values = [
            100 * block.loc[block["role"] == role, column].astype(bool).mean()
            for role in ROLES
        ]
        ax.bar(
            x + (offset - 1) * width,
            values,
            width=width,
            color=LOOP_COL[label],
            edgecolor="none",
        )
    ax.set_xticks(x, [ROLE_LABEL[role] for role in ROLES])
    ax.set_ylim(0, 100)
    if show_ylabel:
        ax.set_ylabel("Loop-associated genes (%)")
    else:
        ax.set_yticklabels([])


def atac_tests(features):
    rows = []
    for trajectory in TRAJECTORIES:
        block = features[features["trajectory"] == trajectory]
        pc = (
            block[block["role"].isin(["P", "C"])]
            .pivot_table(
                index="pair_id", columns="role",
                values="mean_ATAC_TSS_1kb", aggfunc="first"
            )
            .dropna()
        )
        statistic, p_value = wilcoxon(pc["P"], pc["C"])
        rows.append(
            {
                "trajectory": trajectory, "comparison": "P vs C",
                "test": "paired Wilcoxon", "n": len(pc),
                "statistic": statistic, "p_raw": p_value,
            }
        )
        reference = block.loc[block["role"] == "S", "mean_ATAC_TSS_1kb"].dropna()
        for role in ["P", "C"]:
            values = block.loc[
                block["role"] == role, "mean_ATAC_TSS_1kb"
            ].dropna()
            statistic, p_value = mannwhitneyu(
                values, reference, alternative="two-sided"
            )
            rows.append(
                {
                    "trajectory": trajectory, "comparison": f"{role} vs Single",
                    "test": "Mann-Whitney U", "n": len(values),
                    "statistic": statistic, "p_raw": p_value,
                }
            )
    return finish_tests(rows)


def compartment_tests(features):
    rows = []
    for trajectory in TRAJECTORIES:
        block = features[features["trajectory"] == trajectory]
        pc = (
            block[block["role"].isin(["P", "C"])]
            .pivot_table(
                index="pair_id", columns="role",
                values="compartment", aggfunc="first"
            )
            .dropna()
        )
        discordant, p_value = paired_binary_test(pc["P"].eq("A"), pc["C"].eq("A"))
        rows.append(
            {
                "trajectory": trajectory, "comparison": "P vs C",
                "test": "exact McNemar", "statistic": discordant,
                "p_raw": p_value,
            }
        )
        reference = block.loc[block["role"] == "S", "compartment"]
        for role in ["P", "C"]:
            observed = block.loc[block["role"] == role, "compartment"]
            table = [
                [(observed == "A").sum(), (observed == "B").sum()],
                [(reference == "A").sum(), (reference == "B").sum()],
            ]
            odds_ratio, p_value = fisher_exact(table)
            rows.append(
                {
                    "trajectory": trajectory, "comparison": f"{role} vs Single",
                    "test": "Fisher exact", "statistic": odds_ratio,
                    "p_raw": p_value,
                }
            )
    return finish_tests(rows)


def tad_tests(features):
    categories = ["Boundary", "Interior", "Outside"]
    rows = []
    for trajectory in TRAJECTORIES:
        block = features[features["trajectory"] == trajectory]
        pc = (
            block[block["role"].isin(["P", "C"])]
            .pivot_table(
                index="pair_id", columns="role",
                values="TAD_location", aggfunc="first"
            )
            .dropna()
        )
        statistic, df, p_value = paired_multicategory_bowker(
            pc["P"], pc["C"], categories
        )
        rows.append(
            {
                "trajectory": trajectory, "comparison": "P vs C",
                "test": "Bowker symmetry", "statistic": statistic,
                "df": df, "p_raw": p_value,
            }
        )
        reference = block.loc[block["role"] == "S", "TAD_location"]
        for role in ["P", "C"]:
            observed = block.loc[block["role"] == role, "TAD_location"]
            table = np.array(
                [
                    [(observed == category).sum() for category in categories],
                    [(reference == category).sum() for category in categories],
                ]
            )
            statistic, p_value, df, _ = chi2_contingency(table)
            rows.append(
                {
                    "trajectory": trajectory, "comparison": f"{role} vs Single",
                    "test": "Chi-square", "statistic": statistic,
                    "df": df, "p_raw": p_value,
                }
            )
    return finish_tests(rows)


def loop_tests(features):
    rows = []
    for trajectory in TRAJECTORIES:
        block = features[features["trajectory"] == trajectory]
        for label, column in [("Any loop", "any_loop"), ("GGL", "GGL"), ("DGL", "DGL")]:
            pc = (
                block[block["role"].isin(["P", "C"])]
                .pivot_table(
                    index="pair_id", columns="role",
                    values=column, aggfunc="first"
                )
                .dropna()
            )
            discordant, p_value = paired_binary_test(pc["P"], pc["C"])
            rows.append(
                {
                    "trajectory": trajectory, "feature": label,
                    "comparison": "P vs C", "test": "exact McNemar",
                    "statistic": discordant, "p_raw": p_value,
                }
            )
            reference = block.loc[block["role"] == "S", column].astype(bool)
            for role in ["P", "C"]:
                observed = block.loc[block["role"] == role, column].astype(bool)
                table = [
                    [observed.sum(), (~observed).sum()],
                    [reference.sum(), (~reference).sum()],
                ]
                odds_ratio, p_value = fisher_exact(table)
                rows.append(
                    {
                        "trajectory": trajectory, "feature": label,
                        "comparison": f"{role} vs Single", "test": "Fisher exact",
                        "statistic": odds_ratio, "p_raw": p_value,
                    }
                )
    return finish_tests(rows)


def main():
    features = pd.read_csv(SRC / "Figure4_gene_features.tsv.gz", sep="\t")
    profile = pd.read_csv(SRC / "Figure4_ATAC_profile_summary.tsv.gz", sep="\t")
    profile = profile[profile["position"].between(-2500, 2500)].copy()
    for column in ["any_loop", "GGL", "DGL"]:
        features[column] = features[column].astype(str).str.lower().eq("true")

    fig = plt.figure(figsize=(7.2, 8.6))
    grid = fig.add_gridspec(
        4, 4, height_ratios=[1.32, 1, 1, 1],
        left=0.10, right=0.985, top=0.92, bottom=0.065,
        wspace=0.32, hspace=0.42
    )
    for ci, trajectory in enumerate(TRAJECTORIES):
        block = features[features["trajectory"] == trajectory]
        ax_a = fig.add_subplot(grid[0, ci])
        plot_atac(ax_a, profile, trajectory, ci == 0)
        counts = block.groupby("role")["gene"].nunique()
        ax_a.set_title(
            f"{trajectory}\nP/C/Single n={counts.get('P',0)}/{counts.get('C',0)}/{counts.get('S',0)}",
            fontsize=7.5, weight="bold", pad=5, linespacing=1.15
        )
        ax_b = fig.add_subplot(grid[1, ci])
        plot_compartment(ax_b, block, ci == 0)
        ax_c = fig.add_subplot(grid[2, ci])
        plot_tad(ax_c, block, ci == 0)
        ax_d = fig.add_subplot(grid[3, ci])
        plot_loops(ax_d, block, ci == 0)

    fig.text(0.018, 0.958, "A", fontsize=13, weight="bold", va="top")
    for label, row in zip(["B", "C", "D"], [1, 2, 3]):
        fig.text(
            0.018, grid[row, 0].get_position(fig).y1 + 0.018,
            label, fontsize=13, weight="bold", va="top"
        )

    line_handles = [
        Line2D([0], [0], color=RCOL[role], lw=1.6, label=ROLE_LABEL[role]) for role in ROLES
    ]
    fig.legend(
        handles=line_handles, ncol=3, frameon=False,
        loc="upper center", bbox_to_anchor=(0.55, 0.998)
    )
    fig.legend(
        handles=[Patch(facecolor=COMP_COL[x], label=x) for x in ["A", "B"]],
        ncol=2, frameon=False, loc="center right",
        bbox_to_anchor=(0.985, grid[1, 3].get_position(fig).y1 + 0.025)
    )
    fig.legend(
        handles=[Patch(facecolor=TAD_COL[x], label=x) for x in ["Boundary", "Interior", "Outside"]],
        ncol=3, frameon=False, loc="center right",
        bbox_to_anchor=(0.985, grid[2, 3].get_position(fig).y1 + 0.025)
    )
    fig.legend(
        handles=[Patch(facecolor=LOOP_COL[x], label=x) for x in ["Any loop", "GGL", "DGL"]],
        ncol=3, frameon=False, loc="center right",
        bbox_to_anchor=(0.985, grid[3, 3].get_position(fig).y1 + 0.025)
    )

    base = OUT / "Figure4_hapA_two_method_consensus"
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=300)
    fig.savefig(
        base.with_suffix(".tiff"), dpi=600,
        pil_kwargs={"compression": "tiff_lzw"}
    )
    plt.close(fig)

    atac = atac_tests(features)
    comp = compartment_tests(features)
    tad = tad_tests(features)
    loops = loop_tests(features)
    atac.to_csv(OUT / "Figure4A_ATAC_tests_BH.tsv", sep="\t", index=False)
    comp.to_csv(OUT / "Figure4B_compartment_tests_BH.tsv", sep="\t", index=False)
    tad.to_csv(OUT / "Figure4C_TAD_tests_BH.tsv", sep="\t", index=False)
    loops.to_csv(OUT / "Figure4D_loop_tests_BH.tsv", sep="\t", index=False)

    comp_source = (
        features[features["compartment"].isin(["A", "B"])]
        .groupby(["trajectory", "role", "compartment"]).size()
        .rename("n").reset_index()
    )
    tad_source = (
        features.groupby(["trajectory", "role", "TAD_location"]).size()
        .rename("n").reset_index()
    )
    loop_source = (
        features.groupby(["trajectory", "role"])[["any_loop", "GGL", "DGL"]]
        .agg(["sum", "count"]).reset_index()
    )
    comp_source.to_csv(OUT / "Figure4B_compartment_source.tsv", sep="\t", index=False)
    tad_source.to_csv(OUT / "Figure4C_TAD_source.tsv", sep="\t", index=False)
    loop_source.to_csv(OUT / "Figure4D_loop_source.tsv", sep="\t", index=False)
    profile.to_csv(OUT / "Figure4A_ATAC_profile_source.tsv.gz", sep="\t", index=False)

    summary = {
        "trajectory_pairs": 1290,
        "reference_genes": int(features.loc[features["role"] == "S", "gene"].nunique()),
        "ATAC_tests": len(atac),
        "compartment_tests": len(comp),
        "TAD_tests": len(tad),
        "loop_tests": len(loops),
        "reference_definition": (
            "Single = pear gene retained as a one-to-one pear-peach single-copy "
            "ortholog; no peach Hi-C/ATAC data are used."
        ),
    }
    (OUT / "Figure4_analysis_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    image = Image.open(base.with_suffix(".tiff"))
    assert image.info["dpi"][0] >= 599
    assert "<text" in base.with_suffix(".svg").read_text(encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()






