#!/usr/bin/env python3
"""Draw Supplementary Figure S1 after a support-based tail filter."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data" / "supplementary_figures" / "contact_decay_full"
OUT = ROOT / "output" / "supplementary_figure_s1"
OUT.mkdir(parents=True, exist_ok=True)

MIN_POSSIBLE_BIN_PAIRS = 500
SERIES = [
    ("Replicate 1", "contact_decay.rep1.tsv", "#4C78A8"),
    ("Replicate 2", "contact_decay.rep2.tsv", "#E07A5F"),
    ("Merged", "contact_decay.merged.tsv", "#2A9D8F"),
]

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


def load_one(label: str, filename: str) -> pd.DataFrame:
    data = pd.read_csv(SRC / filename, sep="\t")
    required = {"distance_bp", "possible_bin_pairs", "mean_contact_per_possible_pair"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"{filename}: missing columns {sorted(missing)}")
    data = data.copy()
    data.insert(0, "series", label)
    data["display_in_primary_figure"] = (
        data["possible_bin_pairs"] >= MIN_POSSIBLE_BIN_PAIRS
    )
    return data


def main() -> None:
    full = pd.concat(
        [load_one(label, filename) for label, filename, _ in SERIES],
        ignore_index=True,
    )
    shown = full.loc[full["display_in_primary_figure"]].copy()
    full.to_csv(OUT / "Supplementary_Figure_S1_source_full.tsv", sep="\t", index=False)
    shown.to_csv(OUT / "Supplementary_Figure_S1_source.tsv", sep="\t", index=False)

    rows = []
    fig, ax = plt.subplots(figsize=(6.9, 4.2))
    for label, _, color in SERIES:
        block = shown.loc[shown["series"] == label].sort_values("distance_bp")
        ax.plot(
            block["distance_bp"] / 1e6,
            block["mean_contact_per_possible_pair"],
            color=color,
            linewidth=1.45,
            label=label,
        )
        rho, p_value = spearmanr(
            block["distance_bp"], block["mean_contact_per_possible_pair"]
        )
        rows.append(
            {
                "series": label,
                "minimum_possible_bin_pairs": MIN_POSSIBLE_BIN_PAIRS,
                "displayed_distance_bins": len(block),
                "maximum_displayed_distance_bp": int(block["distance_bp"].max()),
                "spearman_rho_distance_vs_mean_contact": rho,
                "spearman_p": p_value,
            }
        )

    ax.set_xlabel("Genomic distance (Mb)")
    ax.set_ylabel("Mean contact per possible bin pair")
    ax.set_yscale("log")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.legend(frameon=False, ncol=3, loc="upper right")
    fig.tight_layout(pad=0.8)
    base = OUT / "Supplementary_Figure_S1"
    for suffix, kwargs in [
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 600}),
        (".tiff", {"dpi": 600}),
    ]:
        fig.savefig(base.with_suffix(suffix), bbox_inches="tight", **kwargs)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(
        OUT / "Supplementary_Figure_S1_summary.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
