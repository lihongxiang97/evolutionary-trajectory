#!/usr/bin/env python3
"""Figure 1 v5: separate cis/trans Circos panels with external legend."""
from pathlib import Path
import gzip
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.65
plt.rcParams["legend.frameon"] = False

import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from make_figure1 import read_bg, add_label, crop_white  # noqa: E402
from make_figure1_v4 import matrix_crop, clean, vertical_gradient  # noqa: E402

SRC, TRACK, OUT = HERE / "source", HERE / "source" / "tracks", HERE / "output"
RED, BLUE = "#C94B44", "#4F75B5"
GOLD, TEAL, PURPLE, DARK = "#D5A11E", "#368C8C", "#8A6BA8", "#3E3E3E"
CHR_COLORS = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#6B8EAD", "#D99058", "#C76C73", "#67A39B", "#7E9D62",
    "#C7A64B", "#9578A6",
]


def chromosome_layout(sizes, gap_fraction=0.004):
    total = sizes.length.sum()
    gap = 2*np.pi*gap_fraction
    usable = 2*np.pi-gap*len(sizes)
    starts, ends, cursor = {}, {}, np.pi/2
    for row in sizes.itertuples():
        span = usable*row.length/total
        starts[row.chrom], ends[row.chrom] = cursor, cursor+span
        cursor += span+gap
    return starts, ends


def theta_for(chrom, pos, size_map, starts, ends):
    return starts[chrom] + (ends[chrom]-starts[chrom])*pos/size_map[chrom]


def line_track(ax, df, size_map, starts, ends, r0, width, color):
    vals = df.value.to_numpy(float)
    lo, hi = np.nanpercentile(vals, [3, 97])
    scale = max(hi-lo, 1e-12)
    for chrom, block in df.groupby("chrom", sort=False):
        if chrom not in starts:
            continue
        mid = (block.start.to_numpy()+block.end.to_numpy())/2
        th = starts[chrom] + (ends[chrom]-starts[chrom])*mid/size_map[chrom]
        v = np.clip((block.value.to_numpy(float)-lo)/scale, 0, 1)
        r = r0+width*v
        ax.plot(r*np.cos(th), r*np.sin(th), color=color, lw=0.42, alpha=0.95)


def pc1_track(ax, df, size_map, starts, ends, r0, width):
    vals = df.value.to_numpy(float)
    scale = max(np.nanpercentile(np.abs(vals), 97), 1e-12)
    seg_pos, seg_neg = [], []
    for row in df.itertuples():
        if row.chrom not in starts:
            continue
        th = theta_for(row.chrom, (row.start+row.end)/2, size_map, starts, ends)
        rr = width*min(abs(row.value)/scale, 1)
        end = r0+rr if row.value >= 0 else r0-rr
        seg = [(r0*np.cos(th), r0*np.sin(th)), (end*np.cos(th), end*np.sin(th))]
        (seg_pos if row.value >= 0 else seg_neg).append(seg)
    ax.add_collection(LineCollection(seg_pos, colors=RED, linewidths=0.35, alpha=0.9))
    ax.add_collection(LineCollection(seg_neg, colors=BLUE, linewidths=0.35, alpha=0.9))


def draw_circos(ax, sizes, mode):
    starts, ends = chromosome_layout(sizes)
    size_map = dict(zip(sizes.chrom, sizes.length))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.12, 1.12)
    ax.set_ylim(-1.12, 1.12)

    # All elements use this one Cartesian coordinate system.
    for i, row in enumerate(sizes.itertuples()):
        th = np.linspace(starts[row.chrom], ends[row.chrom], 140)
        ax.plot(np.cos(th), np.sin(th), lw=7.0, color=CHR_COLORS[i],
                solid_capstyle="butt")
        mid = (starts[row.chrom]+ends[row.chrom])/2
        ax.text(1.075*np.cos(mid), 1.075*np.sin(mid),
                row.chrom.replace("Chr", ""), fontsize=5.0,
                ha="center", va="center", color=DARK)

    # Thin neutral separators between quantitative tracks.
    tt = np.linspace(0, 2*np.pi, 720)
    for r in (0.935, 0.82, 0.70, 0.58):
        ax.plot(r*np.cos(tt), r*np.sin(tt), color="#E3E3E3", lw=0.35)

    pc1 = read_bg(TRACK/"PC1.50kb.bedgraph")
    atac = read_bg(TRACK/"ATAC.mean.50kb.bedgraph")
    genes = read_bg(TRACK/"gene_density.50kb.bedgraph")
    gc = read_bg(TRACK/"GC.50kb.bedgraph")
    pc1_track(ax, pc1, size_map, starts, ends, 0.905, 0.042)
    line_track(ax, atac, size_map, starts, ends, 0.77, 0.050, GOLD)
    line_track(ax, genes, size_map, starts, ends, 0.65, 0.050, TEAL)
    line_track(ax, gc, size_map, starts, ends, 0.53, 0.050, PURPLE)

    with gzip.open(TRACK/"top_loops.for_circos.tsv.gz", "rt") as fh:
        loops = pd.read_csv(fh, sep="\t")
    if mode == "cis":
        loops = loops.loc[loops.chrom1 == loops.chrom2]
        loop_color = RED
    else:
        loops = loops.loc[loops.chrom1 != loops.chrom2]
        loop_color = BLUE

    Path = mpath.Path
    for row in loops.itertuples(index=False):
        a1 = theta_for(row.chrom1, (row.start1+row.end1)/2,
                       size_map, starts, ends)
        a2 = theta_for(row.chrom2, (row.start2+row.end2)/2,
                       size_map, starts, ends)
        p1 = np.array([0.47*np.cos(a1), 0.47*np.sin(a1)])
        p2 = np.array([0.47*np.cos(a2), 0.47*np.sin(a2)])
        path = Path([p1, p1*0.10, p2*0.10, p2],
                    [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(mpatches.PathPatch(
            path, facecolor="none", edgecolor=loop_color,
            lw=0.32 if mode == "cis" else 0.38,
            alpha=0.14 if mode == "cis" else 0.18))


def main():
    a_img = matrix_crop(SRC/"A_whole_genome_50kb_600dpi.png")
    b_img = matrix_crop(SRC/"B_Chr12_10kb_600dpi.png")
    c_img = crop_white(SRC/"C_Chr12_5_10Mb_hicPlotTADs_v4.png", pad=4)
    sizes = pd.read_csv(TRACK/"chrom.sizes", sep="\t", header=None,
                        names=["chrom", "length"])
    sizes["order"] = sizes.chrom.str.replace("Chr", "", regex=False).astype(int)
    sizes = sizes.sort_values("order")

    fig = plt.figure(figsize=(7.2, 9.1), facecolor="white")
    outer = fig.add_gridspec(
        3, 2, height_ratios=[3.0, 2.65, 3.0],
        left=0.075, right=0.97, bottom=0.045, top=0.985,
        wspace=0.24, hspace=0.17)

    # a
    ag = outer[0, 0].subgridspec(1, 2, width_ratios=[1, 0.035], wspace=0.055)
    ax_a = fig.add_subplot(ag[0])
    total = sizes.length.sum()/1e6
    ax_a.imshow(a_img, extent=(0, total, total, 0),
                interpolation="nearest", aspect="equal")
    centers = (sizes.length.cumsum()-sizes.length/2)/1e6
    labels = sizes.chrom.str.replace("Chr", "", regex=False)
    ax_a.set_xticks(centers, labels, rotation=90, fontsize=4.4)
    ax_a.set_yticks(centers, labels, fontsize=4.4)
    ax_a.set_xlabel("Chromosome", labelpad=1)
    ax_a.set_ylabel("Chromosome", labelpad=1)
    clean(ax_a)
    add_label(ax_a, "a", x=-0.12, y=1.01)
    cax = fig.add_subplot(ag[1])
    p = cax.get_position()
    cax.set_position([p.x0, p.y0+p.height*0.18, p.width, p.height*0.64])
    vertical_gradient(cax)

    # b
    bg = outer[0, 1].subgridspec(2, 1, height_ratios=[1, 0.18], hspace=0.035)
    ax_b = fig.add_subplot(bg[0])
    chr12_len = sizes.loc[sizes.chrom == "Chr12", "length"].iloc[0]/1e6
    ax_b.imshow(b_img, extent=(0, chr12_len, chr12_len, 0),
                interpolation="nearest", aspect="equal")
    ax_b.set_ylabel("Chr12 position (Mb)", labelpad=2)
    ax_b.tick_params(labelbottom=False)
    clean(ax_b)
    add_label(ax_b, "b", x=-0.12, y=1.01)
    ax_pc = fig.add_subplot(bg[1], sharex=ax_b)
    pc = read_bg(TRACK/"PC1.50kb.bedgraph")
    pc = pc.loc[pc.chrom == "Chr12"].sort_values("start")
    x = (pc.start+pc.end)/2/1e6
    ax_pc.bar(x, pc.value, width=0.05,
              color=np.where(pc.value >= 0, RED, BLUE), linewidth=0)
    ax_pc.axhline(0, color=DARK, lw=0.4)
    ax_pc.set_xlim(0, chr12_len)
    ax_pc.set_ylabel("PC1", fontsize=6, labelpad=3)
    ax_pc.set_xlabel("Chr12 position (Mb)", labelpad=1)
    clean(ax_pc)

    # c
    ax_c = fig.add_subplot(outer[1, :])
    ax_c.imshow(c_img, interpolation="none", aspect="equal")
    ax_c.axis("off")
    add_label(ax_c, "c", x=-0.025, y=1.015)

    # d/e with an external shared legend.
    dg = outer[2, :].subgridspec(2, 2, height_ratios=[1, 0.12],
                                 hspace=0.01, wspace=0.10)
    ax_d = fig.add_subplot(dg[0, 0])
    draw_circos(ax_d, sizes, "cis")
    add_label(ax_d, "d", x=-0.02, y=0.99)
    ax_d.text(0.5, -0.035, "Intra-chromosomal loops",
              transform=ax_d.transAxes, ha="center", fontsize=6.5)
    ax_e = fig.add_subplot(dg[0, 1])
    draw_circos(ax_e, sizes, "trans")
    add_label(ax_e, "e", x=-0.02, y=0.99)
    ax_e.text(0.5, -0.035, "Inter-chromosomal loops",
              transform=ax_e.transAxes, ha="center", fontsize=6.5)
    lax = fig.add_subplot(dg[1, :])
    lax.axis("off")
    handles = [
        mpatches.Patch(color=RED, label="A compartment (PC1+)"),
        mpatches.Patch(color=BLUE, label="B compartment (PC1−)"),
        Line2D([0], [0], color=GOLD, lw=2, label="ATAC"),
        Line2D([0], [0], color=TEAL, lw=2, label="Gene density"),
        Line2D([0], [0], color=PURPLE, lw=2, label="GC content"),
        Line2D([0], [0], color=RED, lw=1.2, label="cis loops"),
        Line2D([0], [0], color=BLUE, lw=1.2, label="trans loops"),
    ]
    lax.legend(handles=handles, loc="center", ncol=4, fontsize=5.5,
               columnspacing=1.2, handlelength=1.8, handletextpad=0.45)

    # Match B/PC1 physical boundaries exactly.
    fig.canvas.draw()
    bp, pp = ax_b.get_position(), ax_pc.get_position()
    ax_pc.set_position([bp.x0, pp.y0, bp.width, pp.height])

    for fmt in ("svg", "pdf", "png", "tiff"):
        kw = {"bbox_inches": "tight", "facecolor": "white"}
        if fmt == "png":
            kw["dpi"] = 600
        elif fmt == "tiff":
            kw.update(dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        fig.savefig(OUT/f"Figure1_hapA.{fmt}", **kw)
    plt.close(fig)


if __name__ == "__main__":
    main()
