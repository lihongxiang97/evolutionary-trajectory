#!/usr/bin/env python3
"""Rebuild manuscript Figure 1 from the hapA 3D-genome reanalysis."""
from pathlib import Path
import gzip
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["legend.frameon"] = False

import matplotlib.path as mpath
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT / "source_data" / "main_figures" / "figure1"
OUT = ROOT / "output" / "figure1"
TRACK = SRC / "tracks"
OUT.mkdir(parents=True, exist_ok=True)

RED = "#C94B44"
BLUE = "#4F75B5"
GOLD = "#D5A11E"
TEAL = "#368C8C"
PURPLE = "#8A6BA8"
DARK = "#3E3E3E"
LIGHT = "#E5E5E5"


def read_bg(path):
    return pd.read_csv(path, sep="\t", header=None,
                       names=["chrom", "start", "end", "value"])


def crop_white(path, pad=8):
    im = Image.open(path).convert("RGB")
    a = np.asarray(im)
    mask = np.any(a < 248, axis=2)
    yy, xx = np.where(mask)
    if len(xx):
        box = (max(0, xx.min()-pad), max(0, yy.min()-pad),
               min(im.width, xx.max()+pad), min(im.height, yy.max()+pad))
        im = im.crop(box)
    return np.asarray(im)


def add_label(ax, text, x=-0.06, y=1.03):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left")


def chromosome_layout(sizes, gap_fraction=0.004):
    total = sizes["length"].sum()
    gap = 2 * np.pi * gap_fraction
    usable = 2 * np.pi - gap * len(sizes)
    starts, ends, cursor = {}, {}, np.pi / 2
    for row in sizes.itertuples():
        span = usable * row.length / total
        starts[row.chrom], ends[row.chrom] = cursor, cursor + span
        cursor += span + gap
    return starts, ends


def angle(chrom, pos, sizes, starts, ends):
    length = sizes.loc[sizes.chrom == chrom, "length"].iloc[0]
    return starts[chrom] + (ends[chrom] - starts[chrom]) * pos / length


def ring_track(ax, df, sizes, starts, ends, r0, width, color,
               signed=False, linewidth=0.45):
    vals = df["value"].to_numpy(float)
    scale = np.nanpercentile(np.abs(vals), 97)
    scale = scale if scale > 0 else 1.0
    for chrom, block in df.groupby("chrom", sort=False):
        if chrom not in starts:
            continue
        theta = np.array([angle(chrom, (s+e)/2, sizes, starts, ends)
                          for s, e in zip(block.start, block.end)])
        v = np.clip(block.value.to_numpy(float) / scale, -1, 1)
        if signed:
            colors = np.where(v >= 0, RED, BLUE)
            bottom = np.where(v >= 0, r0, r0 - width * np.abs(v))
            height = width * np.abs(v)
            ax.bar(theta, height, width=(ends[chrom]-starts[chrom])/len(block)*1.02,
                   bottom=bottom, color=colors, linewidth=0)
        else:
            rr = r0 + width * np.clip(v, 0, 1)
            ax.plot(theta, rr, color=color, lw=linewidth, alpha=0.9)
            ax.fill_between(theta, r0, rr, color=color, alpha=0.35, linewidth=0)


def draw_circos(ax):
    sizes = pd.read_csv(TRACK / "chrom.sizes", sep="\t", header=None,
                        names=["chrom", "length"])
    sizes["order"] = sizes.chrom.str.replace("Chr", "", regex=False).astype(int)
    sizes = sizes.sort_values("order")
    starts, ends = chromosome_layout(sizes)
    ax.set_aspect("equal")
    ax.axis("off")
    # Chromosome ideograms and labels.
    for i, row in enumerate(sizes.itertuples()):
        t = np.linspace(starts[row.chrom], ends[row.chrom], 100)
        ax.plot(np.cos(t), np.sin(t), lw=7,
                color="#A7A7A7" if i % 2 else "#555555",
                solid_capstyle="butt")
        mid = (starts[row.chrom] + ends[row.chrom]) / 2
        ax.text(1.08*np.cos(mid), 1.08*np.sin(mid), row.chrom.replace("Chr", ""),
                fontsize=5.8, ha="center", va="center")

    # Use a hidden polar inset for quantitative rings.
    polar = ax.inset_axes([0, 0, 1, 1], projection="polar")
    polar.patch.set_alpha(0)
    polar.set_theta_direction(1)
    polar.set_theta_offset(0)
    polar.set_ylim(0, 1.05)
    polar.axis("off")
    pc1 = read_bg(TRACK / "PC1.50kb.bedgraph")
    atac = read_bg(TRACK / "ATAC.mean.50kb.bedgraph")
    genes = read_bg(TRACK / "gene_density.50kb.bedgraph")
    gc = read_bg(TRACK / "GC.50kb.bedgraph")
    ring_track(polar, pc1, sizes, starts, ends, 0.91, 0.055, RED, signed=True)
    ring_track(polar, atac, sizes, starts, ends, 0.79, 0.060, GOLD)
    ring_track(polar, genes, sizes, starts, ends, 0.68, 0.060, TEAL)
    gc["value"] = (gc["value"] - 0.5).clip(lower=0)
    ring_track(polar, gc, sizes, starts, ends, 0.57, 0.060, PURPLE)

    # Loop chords in Cartesian coordinates.
    Path = mpath.Path
    with gzip.open(TRACK / "top_loops.for_circos.tsv.gz", "rt") as fh:
        loops = pd.read_csv(fh, sep="\t")
    for row in loops.itertuples(index=False):
        m1, m2 = (row.start1 + row.end1)/2, (row.start2 + row.end2)/2
        a1 = angle(row.chrom1, m1, sizes, starts, ends)
        a2 = angle(row.chrom2, m2, sizes, starts, ends)
        p1 = np.array([0.52*np.cos(a1), 0.52*np.sin(a1)])
        p2 = np.array([0.52*np.cos(a2), 0.52*np.sin(a2)])
        # Pull control points toward the centre to obtain Circos-style chords.
        c1, c2 = p1 * 0.12, p2 * 0.12
        path = Path([p1, c1, c2, p2],
                    [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        col = RED if row.chrom1 == row.chrom2 else BLUE
        ax.add_patch(mpatches.PathPatch(path, facecolor="none", edgecolor=col,
                                       lw=0.28, alpha=0.12))
    ax.text(-1.34, 0.68, "PC1", color=DARK, fontsize=6)
    ax.text(-1.34, 0.56, "ATAC", color=GOLD, fontsize=6)
    ax.text(-1.34, 0.44, "Gene density", color=TEAL, fontsize=6)
    ax.text(-1.34, 0.32, "GC", color=PURPLE, fontsize=6)
    ax.text(-1.34, 0.20, "cis loops", color=RED, fontsize=6)
    ax.text(-1.34, 0.08, "trans loops", color=BLUE, fontsize=6)
    ax.set_xlim(-1.42, 1.18)
    ax.set_ylim(-1.15, 1.15)


def main():
    required = [SRC / "A_whole_genome_50kb.png", SRC / "B_Chr12_10kb.png",
                SRC / "C_Chr12_5_10Mb_TAD.png"]
    missing = [str(x) for x in required if not x.exists() or x.stat().st_size == 0]
    if missing:
        raise SystemExit("Missing heatmap panels: " + ", ".join(missing))

    fig = plt.figure(figsize=(7.2, 9.1), facecolor="white")
    outer = fig.add_gridspec(3, 2, height_ratios=[2.55, 2.10, 4.0],
                             width_ratios=[1, 1], hspace=0.26, wspace=0.14)

    ax_a = fig.add_subplot(outer[0, 0])
    ax_a.imshow(crop_white(required[0]))
    ax_a.axis("off")
    ax_a.set_title("Genome-wide contacts (50 kb)", fontsize=8, pad=2)
    add_label(ax_a, "a")

    bgs = outer[0, 1].subgridspec(2, 1, height_ratios=[5, 1.05], hspace=0.02)
    ax_b = fig.add_subplot(bgs[0])
    ax_b.imshow(crop_white(required[1]))
    ax_b.axis("off")
    ax_b.set_title("Chr12 contacts and compartments (10 kb / 50 kb)", fontsize=8, pad=2)
    add_label(ax_b, "b")
    ax_pc = fig.add_subplot(bgs[1])
    pc = read_bg(TRACK / "PC1.50kb.bedgraph")
    pc = pc.loc[pc.chrom == "Chr12"].sort_values("start")
    x = (pc.start + pc.end) / 2 / 1e6
    colors = np.where(pc.value >= 0, RED, BLUE)
    ax_pc.bar(x, pc.value, width=0.05, color=colors, linewidth=0)
    ax_pc.axhline(0, color=DARK, lw=0.4)
    ax_pc.set_xlim(0, pc.end.max()/1e6)
    ax_pc.set_ylabel("PC1", fontsize=6)
    ax_pc.set_xlabel("Chr12 position (Mb)", fontsize=6)
    ax_pc.tick_params(labelsize=5, length=2)
    ax_pc.spines[["top", "right"]].set_visible(False)

    cgs = outer[1, :].subgridspec(2, 1, height_ratios=[4.8, 1.1], hspace=0.03)
    ax_c = fig.add_subplot(cgs[0])
    ax_c.imshow(crop_white(required[2]))
    ax_c.axis("off")
    ax_c.set_title("TAD organization in Chr12:5–10 Mb (10 kb)", fontsize=8, pad=2)
    add_label(ax_c, "c", x=-0.025)
    ax_at = fig.add_subplot(cgs[1])
    at = read_bg(TRACK / "ATAC.mean.Chr12.5_10Mb.10kb.bedgraph")
    xx = (at.start + at.end) / 2 / 1e6
    ax_at.fill_between(xx, 0, at.value, color=GOLD, alpha=0.75, lw=0)
    ax_at.plot(xx, at.value, color="#A97900", lw=0.45)
    ax_at.set_xlim(5, 10)
    ax_at.set_ylabel("ATAC\n(RPGC)", fontsize=6)
    ax_at.set_xlabel("Chr12 position (Mb)", fontsize=6)
    ax_at.tick_params(labelsize=5, length=2)
    ax_at.spines[["top", "right"]].set_visible(False)

    ax_d = fig.add_subplot(outer[2, :])
    draw_circos(ax_d)
    ax_d.set_title("Genome-wide chromatin loops and genomic features", fontsize=8, pad=0)
    add_label(ax_d, "d", x=-0.01, y=0.97)

    fig.suptitle("Three-dimensional organization of the pear hapA genome",
                 fontsize=10, fontweight="bold", y=0.995)
    for fmt in ("svg", "pdf", "png", "tiff"):
        kwargs = {"bbox_inches": "tight"}
        if fmt in ("png", "tiff"):
            kwargs["dpi"] = 600 if fmt == "tiff" else 300
        if fmt == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(OUT / f"Figure1_hapA.{fmt}", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
