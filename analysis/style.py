"""Shared figure style for the README analysis (spec #34, ticket #39): sizes in millimetres, fonts of at least 7 pt,
the Okabe-Ito colour-blind-safe palette, editable text, deterministic files. Every figure script imports this module."""
import matplotlib
import matplotlib.pyplot as plt

OKABE_ITO = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
# one colour per database, the same in every figure (colour is never the only cue: names are labelled too)
DATABASE_COLOURS = {"card": "#009E73", "ncbi": "#0072B2", "resfinder": "#D55E00"}
FULL = "#D55E00"      # full-name regression
ALLELE_FREE = "#E69F00"  # allele-free regression
GREY = "#555555"
FONT_PT = 7


def apply_style():
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": FONT_PT, "axes.labelsize": FONT_PT, "axes.titlesize": FONT_PT, "xtick.labelsize": FONT_PT,
        "ytick.labelsize": FONT_PT, "legend.fontsize": FONT_PT, "figure.titlesize": FONT_PT + 1,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "lines.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,   # editable text
        "svg.fonttype": "none", "svg.hashsalt": "chamredb",  # editable text and reproducible ids
    })


def new_figure(width_mm, height_mm):
    """A figure of the final size in millimetres (design at final size)."""
    apply_style()
    return plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))


def axes_mm(fig, left, bottom, width, height):
    """Axes placed in millimetres from the bottom left corner of the figure."""
    w, h = fig.get_size_inches() * 25.4
    return fig.add_axes([left / w, bottom / h, width / w, height / h])


def save_figure(fig, base_path):
    """Write base_path.svg and base_path.png (300 dpi) with no date or software metadata (byte-identical re-runs)."""
    fig.savefig(f"{base_path}.svg", metadata={"Date": None})
    fig.savefig(f"{base_path}.png", dpi=300, metadata={"Software": None})
