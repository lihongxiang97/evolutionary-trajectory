#!/usr/bin/env python3
"""Figure 1 revision: vertical legend, aligned tracks, TAD lines, fixed Circos."""
from pathlib import Path
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

import numpy as np
import pandas as pd
from PIL import Image
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from make_figure1_base import read_bg, draw_circos, add_label, crop_white  # noqa: E402

SRC = ROOT / "source_data" / "main_figures" / "figure1"
TRACK = SRC / "tracks"
OUT = ROOT / "output" / "figure1"
OUT.mkdir(parents=True, exist_ok=True)
RED, BLUE, DARK = "#C94B44", "#4F75B5", "#3E3E3E"


def matrix_crop(path):
    im = Image.open(path).convert("RGB")
    w, h = im.size
    box = (round(w * 452/2400), round(h * 272/2100),
           round(w * 1886/2400), round(h * 1709/2100))
    return np.asarray(im.crop(box))


def clean(ax):
    ax.tick_params(length=2, width=0.5, labelsize=5)
    ax.spines[["top", "right"]].set_visible(False)


def vertical_gradient(ax):
    cmap = LinearSegmentedColormap.from_list(
        "contact", ["#fff7f3", "#fcbba1", "#ef3b2c", "#7f0000"])
    ax.imshow(np.linspace(0, 1, 256)[:, None], aspect="auto",
              cmap=cmap, origin="lower")
    ax.set_xticks([])
    ax.set_yticks([0, 255], ["Low", "High"], fontsize=5)
    ax.tick_params(length=0, pad=2)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("Normalized contact frequency", rotation=270,
                  va="bottom", labelpad=7, fontsize=5.5)
    for spine in ax.spines.values():
        spine.set_visible(False)


def fixed_circos(ax):
    draw_circos(ax)
    # Parent Cartesian circle and inset polar tracks now use identical
    # centre/radius mapping: both span -1.05..1.05 or r=0..1.05.
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    track_names = ["PC1", "ATAC", "Gene density", "GC", "cis loops", "trans loops"]
    y_positions = [0.31, 0.265, 0.22, 0.175, 0.13, 0.085]
    for txt in ax.texts:
        if txt.get_text() in track_names:
            i = track_names.index(txt.get_text())
            txt.set_transform(ax.transAxes)
            txt.set_position((0.015, y_positions[i]))
            txt.set_ha("left")
            txt.set_va("center")
            txt.set_clip_on(False)


def main():
    a_img = matrix_crop(SRC / "A_whole_genome_50kb_600dpi.png")
    b_img = matrix_crop(SRC / "B_Chr12_10kb_600dpi.png")
    c_img = crop_white(SRC / "C_Chr12_5_10Mb_hicPlotTADs_v4.png", pad=4)
    sizes = pd.read_csv(TRACK / "chrom.sizes", sep="\t", header=None,
                        names=["chrom", "length"])
    sizes["order"] = sizes.chrom.str.replace("Chr", "", regex=False).astype(int)
    sizes = sizes.sort_values("order")

    fig = plt.figure(figsize=(7.2, 9.0), facecolor="white")
    outer = fig.add_gridspec(
        3, 2, height_ratios=[3.0, 2.25, 3.25],
        left=0.075, right=0.97, bottom=0.035, top=0.985,
        wspace=0.24, hspace=0.18)

    # a: vertical legend immediately to the right of the matrix.
    ag = outer[0, 0].subgridspec(1, 2, width_ratios=[1, 0.035], wspace=0.055)
    ax_a = fig.add_subplot(ag[0])
    total = sizes.length.sum()/1e6
    ax_a.imshow(a_img, extent=(0, total, total, 0),
                interpolation="nearest", aspect="equal")
    centers = (sizes.length.cumsum() - sizes.length/2)/1e6
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

    # b: merged Chr12 heatmap and PC1; actual axes positions are synchronized below.
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
    pc = read_bg(TRACK / "PC1.50kb.bedgraph")
    pc = pc.loc[pc.chrom == "Chr12"].sort_values("start")
    x = (pc.start + pc.end)/2/1e6
    ax_pc.bar(x, pc.value, width=0.05,
              color=np.where(pc.value >= 0, RED, BLUE), linewidth=0)
    ax_pc.axhline(0, color=DARK, lw=0.4)
    ax_pc.set_xlim(0, chr12_len)
    ax_pc.set_ylabel("PC1", fontsize=6, labelpad=3)
    ax_pc.set_xlabel("Chr12 position (Mb)", labelpad=1)
    clean(ax_pc)

    # c: hicPlotTADs with triangular TAD boundary lines, TAD score and mean ATAC.
    ax_c = fig.add_subplot(outer[1, :])
    ax_c.imshow(c_img, interpolation="nearest", aspect="auto")
    ax_c.axis("off")
    add_label(ax_c, "c", x=-0.025, y=1.015)

    # d: corrected single-centre circular coordinate system.
    ax_d = fig.add_subplot(outer[2, :])
    fixed_circos(ax_d)
    add_label(ax_d, "d", x=0.13, y=0.99)

    # Equalize the physical left/right boundaries of the B heatmap and PC1 axes.
    fig.canvas.draw()
    bp, pp = ax_b.get_position(), ax_pc.get_position()
    ax_pc.set_position([bp.x0, pp.y0, bp.width, pp.height])

    for fmt in ("svg", "pdf", "png", "tiff"):
        kw = {"bbox_inches": "tight", "facecolor": "white"}
        if fmt == "png":
            kw["dpi"] = 300
        elif fmt == "tiff":
            kw.update(dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        fig.savefig(OUT / f"Figure1_hapA.{fmt}", **kw)
    plt.close(fig)


if __name__ == "__main__":
    main()
