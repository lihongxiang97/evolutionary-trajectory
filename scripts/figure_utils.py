from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image


YOUNG = "#E07A5F"
OLD = "#4C78A8"
AGE_COLORS = {"young": YOUNG, "old": OLD}
AGE_LABELS = {"young": "Young", "old": "Old"}


def configure():
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.2,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def significance(p):
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def panel_label(fig, ax, label, x_offset=0.035, y_offset=0.015):
    box = ax.get_position(fig)
    fig.text(
        box.x0 - x_offset,
        box.y1 + y_offset,
        label,
        fontsize=12,
        weight="bold",
        ha="left",
        va="top",
    )


def save_bundle(fig, base, dpi_png=300):
    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=dpi_png, bbox_inches="tight")
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    image = Image.open(base.with_suffix(".tiff"))
    assert image.info["dpi"][0] >= 599
    svg = base.with_suffix(".svg").read_text(encoding="utf-8")
    assert "<text" in svg

