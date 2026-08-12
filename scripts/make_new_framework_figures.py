#!/usr/bin/env python3
"""Create Figure 8 and Supplementary Figures S16-S18 with Python only."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns



COLORS = {
    "neither": "#D9D9D9", "both": "#7A6FAC", "child_only": "#D95F02", "parent_only": "#1B9E77",
    "TD": "#4C78A8", "PD": "#F2A541", "TRD": "#A05195", "primary": "#355C7D", "strict": "#C06C84",
}
STATE_ORDER = ["neither", "both", "child_only", "parent_only"]
STATE_LABEL = {"neither": "Neither", "both": "Both", "child_only": "Child only", "parent_only": "Parent only"}
SET_LABEL = {"all_primary_pairs": "All", "TRD": "TRD", "TRD_peach_polarized": "TRD\npeach-polarized",
             "TRD_multi_outgroup_polarized": "TRD\nmulti-outgroup"}


def style():
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7.2, "axes.labelsize": 7.5, "axes.titlesize": 8,
        "xtick.labelsize": 6.7, "ytick.labelsize": 6.7, "legend.fontsize": 6.5,
        "axes.linewidth": 0.7, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
    })
    sns.set_style("ticks")


def panel_label(ax, label):
    ax.text(-0.22, 1.13, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")



def wilson_ci(successes, total, z=1.959963984540054):
    if total <= 0:
        return (np.nan, np.nan)
    p = successes / total
    den = 1 + z * z / total
    center = (p + z * z / (2 * total)) / den
    half = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / den
    return center - half, center + half
def significance(q):
    if pd.isna(q): return ""
    if q < 0.001: return "q<0.001"
    if q < 0.01: return f"q={q:.3f}"
    return f"q={q:.2f}"


def save_all(fig, outdir, stem):
    for ext in ["pdf", "svg", "png", "tiff"]:
        dpi = 600 if ext == "tiff" else 300
        fig.savefig(outdir / f"{stem}.{ext}", dpi=dpi, bbox_inches="tight")


def main_figure(data_dir, outdir, source_dir):
    counts = pd.read_csv(data_dir / "regulatory_capture_state_counts.tsv", sep="\t")
    direction = pd.read_csv(data_dir / "regulatory_capture_direction_tests.tsv", sep="\t")
    models = pd.read_csv(data_dir / "regulatory_capture_adjusted_models.tsv", sep="\t")
    tad = pd.read_csv(data_dir / "TAD_coresidence_distance_matched.tsv", sep="\t")
    coex = pd.read_csv(data_dir / "TAD_coresidence_coexpression_models.tsv", sep="\t")
    hap = pd.read_csv(data_dir / "hapA_hapB_pair_validation_summary.tsv", sep="\t")

    fig, axes = plt.subplots(2, 3, figsize=(7.09, 5.25), gridspec_kw={"wspace": 0.64, "hspace": 0.82})

    ax = axes[0, 0]
    a = counts[(counts.loop_set == "primary") & counts.analysis_set.isin(["all_primary_pairs", "TRD", "TRD_multi_outgroup_polarized"])].copy()
    order = ["all_primary_pairs", "TRD", "TRD_multi_outgroup_polarized"]
    bottom = np.zeros(len(order))
    for state in STATE_ORDER:
        vals = [a.loc[(a.analysis_set == x) & (a.capture_state == state), "fraction"].iloc[0] * 100 for x in order]
        ax.bar(range(len(order)), vals, bottom=bottom, width=0.72, color=COLORS[state], label=STATE_LABEL[state])
        bottom += vals
    ax.set_xticks(range(len(order)), [SET_LABEL[x] for x in order])
    ax.set_ylabel("Duplicate pairs (%)")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.31), columnspacing=0.9, handlelength=1.2)
    panel_label(ax, "A")

    ax = axes[0, 1]
    b = direction[direction.loop_set == "primary"].copy()
    b = b.set_index("analysis_set").loc[order].reset_index()
    y = np.arange(len(b))
    frac, lo, hi = [], [], []
    for row in b.itertuples():
        value = row.child_fraction_among_discordant
        ci = wilson_ci(row.child_only_n, row.discordant_n)
        frac.append(value); lo.append(ci[0]); hi.append(ci[1])
    ax.errorbar(frac, y, xerr=[np.array(frac)-np.array(lo), np.array(hi)-np.array(frac)], fmt="o",
                color="#333333", ecolor="#777777", capsize=2.5, ms=4)
    ax.axvline(0.5, color="#999999", lw=0.8, ls="--")
    ax.set_yticks(y, [SET_LABEL[x] for x in order])
    ax.set_xlim(0.25, 0.75)
    ax.set_xlabel("Child-only fraction among\ndiscordant pairs")
    ax.invert_yaxis()
    for yi, row in enumerate(b.itertuples()):
        ax.text(0.74, yi, f"n={row.discordant_n}", ha="right", va="center", fontsize=6.2)
    panel_label(ax, "B")

    ax = axes[0, 2]
    c = models[models.term.str.contains("capture_delta", na=False) & models.analysis_set.isin(["all", "TRD", "TRD_multi_outgroup_polarized"])].copy()
    c["label"] = c.analysis_set.map({"all": "All", "TRD": "TRD", "TRD_multi_outgroup_polarized": "TRD multi-outgroup"})
    c["order"] = c["label"].map({"All": 0, "TRD": 1, "TRD multi-outgroup": 2}) + c.loop_set.map({"primary": -0.12, "strict": 0.12})
    for loop_set, sub in c.groupby("loop_set"):
        ax.errorbar(sub.estimate, sub.order, xerr=[sub.estimate-sub.ci_low, sub.ci_high-sub.estimate], fmt="o",
                    color=COLORS[loop_set], ecolor=COLORS[loop_set], capsize=2.5, ms=4, label=loop_set.capitalize())
    ax.axvline(0, color="#777777", lw=0.8, ls="--")
    ax.set_yticks([0, 1, 2], ["All", "TRD", "TRD\nmulti-outgroup"])
    ax.set_xlabel("Adjusted child–parent expression\ncoefficient (95% CI)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.27), ncol=2)
    panel_label(ax, "C")

    ax = axes[1, 0]
    d = tad.set_index("duplication_type").loc[["TD", "PD", "TRD"]].reset_index()
    x = np.arange(3); width = 0.34
    ax.bar(x-width/2, d.observed_same_TAD_fraction*100, width, color=[COLORS[x] for x in d.duplication_type], label="Observed")
    err = np.vstack([(d.matched_random_mean-d.matched_random_ci_low)*100,
                     (d.matched_random_ci_high-d.matched_random_mean)*100])
    ax.bar(x+width/2, d.matched_random_mean*100, width, color="#CFCFCF", label="Distance matched",
           yerr=err, capsize=2, error_kw={"lw": 0.7})
    ax.set_xticks(x, d.duplication_type)
    ax.set_ylabel("Pairs in the same TAD (%)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.25), ncol=2)
    for xi, row in enumerate(d.itertuples()):
        ax.text(xi, max(row.observed_same_TAD_fraction, row.matched_random_mean)*100 + 4,
                significance(row.empirical_p_adj_BH), ha="center", fontsize=6.2)
    ax.set_ylim(0, 108)
    panel_label(ax, "D")

    ax = axes[1, 1]
    e = coex[(coex.term == "same_TAD_int") & coex.analysis_set.isin(["all", "TD", "PD", "TRD"])].copy()
    e = e.set_index("analysis_set").loc[["all", "TD", "PD", "TRD"]].reset_index()
    y = np.arange(len(e))
    ax.errorbar(e.estimate, y, xerr=[e.estimate-e.ci_low, e.ci_high-e.estimate], fmt="o", color="#333333",
                ecolor="#777777", capsize=2.5, ms=4)
    ax.axvline(0, color="#777777", lw=0.8, ls="--")
    ax.set_yticks(y, ["All", "TD", "PD", "TRD"])
    ax.invert_yaxis()
    ax.set_xlabel("Adjusted same-TAD effect on\ncoexpression (Fisher z, 95% CI)")
    panel_label(ax, "E")

    ax = axes[1, 2]
    metric_order = ["capture_delta_exact_agreement", "capture_delta_exact_agreement_strict", "same_TAD_agreement"]
    f = hap.set_index("metric").loc[metric_order].reset_index()
    labels = ["Primary\ncapture", "Strict\ncapture", "Same-TAD\nstatus"]
    bars = ax.bar(range(3), f.agreement_fraction*100, color=[COLORS["primary"], COLORS["strict"], "#4C78A8"], width=0.65)
    ax.set_xticks(range(3), labels)
    ax.set_ylim(0, 105)
    ax.set_ylabel("hapA–hapB agreement (%)")
    for bar, row in zip(bars, f.itertuples()):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+2, f"{bar.get_height():.1f}%\n(n={row.n_pairs})",
                ha="center", va="bottom", fontsize=6.1)
    panel_label(ax, "F")

    for ax in axes.flat:
        sns.despine(ax=ax)
    fig.subplots_adjust(left=0.11, right=0.985, top=0.88, bottom=0.13)
    save_all(fig, outdir, "Figure8_3D_duplicate_mechanisms")

    a.to_csv(source_dir / "Figure8A_capture_states.tsv", sep="\t", index=False)
    b.to_csv(source_dir / "Figure8B_capture_direction.tsv", sep="\t", index=False)
    c.to_csv(source_dir / "Figure8C_capture_models.tsv", sep="\t", index=False)
    d.to_csv(source_dir / "Figure8D_TAD_coresidence.tsv", sep="\t", index=False)
    e.to_csv(source_dir / "Figure8E_coexpression_models.tsv", sep="\t", index=False)
    f.to_csv(source_dir / "Figure8F_cross_haplotype_robustness.tsv", sep="\t", index=False)
    plt.close(fig)


def supplementary(data_dir, outdir, source_dir):
    pair = pd.read_csv(data_dir / "duplicate_pair_3D_features.tsv.gz", sep="\t")
    rows = []
    for loop_set, state_col in [("Primary", "capture_state"), ("Strict", "capture_state_strict")]:
        tmp = pair.groupby(["duplication_type", "trajectory", state_col]).size().rename("n").reset_index()
        tmp = tmp.rename(columns={state_col: "capture_state"})
        tmp["loop_set"] = loop_set
        tmp["fraction"] = tmp["n"] / tmp.groupby(["duplication_type", "trajectory"])["n"].transform("sum")
        rows.append(tmp)
    s16 = pd.concat(rows, ignore_index=True)
    fig, axes = plt.subplots(2, 3, figsize=(7.09, 4.7), sharey=True, gridspec_kw={"hspace": 0.42, "wspace": 0.20})
    trajectories = ["Con", "Neo(C)", "Neo(P)", "Spe"]
    for i, loop_set in enumerate(["Primary", "Strict"]):
        for j, dtype in enumerate(["TD", "PD", "TRD"]):
            ax = axes[i, j]
            sub = s16[(s16.loop_set == loop_set) & (s16.duplication_type == dtype)]
            bottom = np.zeros(4)
            for state in STATE_ORDER:
                vals = []
                for tr in trajectories:
                    hit = sub[(sub.trajectory == tr) & (sub.capture_state == state)]
                    vals.append(hit.fraction.iloc[0]*100 if len(hit) else 0)
                ax.bar(range(4), vals, bottom=bottom, color=COLORS[state], width=0.72,
                       label=STATE_LABEL[state] if i == 0 and j == 0 else None)
                bottom += vals
            ax.set_xticks(range(4), trajectories, rotation=25, ha="right")
            ax.set_title(dtype)
            if j == 0: ax.set_ylabel(f"{loop_set}\nPairs (%)")
            if i == 0 and j == 0:
                ax.legend(frameon=False, ncol=4, loc="upper left", bbox_to_anchor=(0, 1.35), columnspacing=0.8)
            sns.despine(ax=ax)
    save_all(fig, outdir, "Supplementary_Figure_S16_capture_states")
    s16.to_csv(source_dir / "Supplementary_Figure_S16.tsv", sep="\t", index=False)
    plt.close(fig)

    age = pd.read_csv(data_dir / "TAD_coresidence_age_models.tsv", sep="\t")
    desc = age[age.analysis == "descriptive"].copy()
    adj = age[(age.analysis == "adjusted_logistic") & (age.term == "young")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.8), gridspec_kw={"wspace": 0.42})
    ax = axes[0]
    dtypes = ["TD", "PD", "TRD"]; x = np.arange(3); width = 0.34
    for k, age_group in enumerate(["young", "old"]):
        vals = [desc[(desc.duplication_type == d) & (desc.term == age_group)].estimate.iloc[0]*100 for d in dtypes]
        ax.bar(x + (k-0.5)*width, vals, width, color=["#E76F51", "#457B9D"][k], label=age_group.capitalize())
    ax.set_xticks(x, dtypes); ax.set_ylabel("Pairs in the same TAD (%)"); ax.legend(frameon=False)
    panel_label(ax, "A")
    ax = axes[1]
    if len(adj):
        adj = adj.set_index("duplication_type").loc[[x for x in ["TD", "PD"] if x in adj.duplication_type.values]].reset_index()
        y = np.arange(len(adj))
        ax.errorbar(adj.estimate, y, xerr=[adj.estimate-adj.ci_low, adj.ci_high-adj.estimate], fmt="o",
                    color="#333333", ecolor="#777777", capsize=2.5)
        ax.set_yticks(y, adj.duplication_type); ax.invert_yaxis()
    ax.axvline(1, color="#777777", lw=0.8, ls="--"); ax.set_xscale("log")
    ax.set_xlabel("Adjusted same-TAD odds ratio\n(Young vs Old, 95% CI)")
    panel_label(ax, "B")
    for ax in axes: sns.despine(ax=ax)
    save_all(fig, outdir, "Supplementary_Figure_S17_TAD_age")
    age.to_csv(source_dir / "Supplementary_Figure_S17.tsv", sep="\t", index=False)
    plt.close(fig)

    hv = pd.read_csv(data_dir / "hapA_hapB_pair_validation.tsv.gz", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.45), gridspec_kw={"wspace": 0.45})
    configs = [("hapA_capture_delta", "hapB_capture_delta", "Capture delta\n(primary)"),
               ("hapA_capture_delta_strict", "hapB_capture_delta_strict", "Capture delta\n(strict)"),
               ("hapA_same_TAD", "hapB_same_TAD", "Same-TAD status")]
    mats = []
    for ax, (ac, bc, title) in zip(axes, configs):
        labels = [-1, 0, 1] if "delta" in ac else [0, 1]
        tab = pd.crosstab(hv[ac], hv[bc]).reindex(index=labels, columns=labels, fill_value=0)
        sns.heatmap(tab, annot=True, fmt="d", cmap="Blues", cbar=False, square=True, ax=ax,
                    linewidths=0.5, linecolor="white")
        ax.set_xlabel("hapB"); ax.set_ylabel("hapA"); ax.set_title(title)
        mats.append(tab.stack().rename("n").reset_index().assign(metric=title.replace("\n", " ")))
    save_all(fig, outdir, "Supplementary_Figure_S18_haplotype_agreement")
    pd.concat(mats, ignore_index=True).to_csv(source_dir / "Supplementary_Figure_S18.tsv", sep="\t", index=False)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    style()
    data_dir = Path(args.data_dir); outdir = Path(args.outdir)
    source_dir = outdir / "source_data"
    outdir.mkdir(parents=True, exist_ok=True); source_dir.mkdir(parents=True, exist_ok=True)
    main_figure(data_dir, outdir, source_dir)
    supplementary(data_dir, outdir, source_dir)


if __name__ == "__main__":
    main()
