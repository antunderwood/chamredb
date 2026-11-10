"""Figures of the README analysis (spec #34). Display transforms only (binning, a fitted line and its 95% interval);
every statistic comes from the result of analysis.matches.assess and analysis.run.statistics."""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter

from analysis import style

NAME = "figure_name_similarity_vs_identity"
WIDTH_MM, HEIGHT_MM = 180, 150
BINS_X, BINS_Y = 40, 20  # name similarity takes discrete values, so coarser bins in y avoid gappy stripes


def _histogram(values, bins=20):
    counts, edges = np.histogram(values, bins=bins, range=(0, 1))
    return pd.DataFrame({"bin_low": edges[:-1], "bin_high": edges[1:], "matches": counts})


def _fit_band(x, y, fit, grid):
    """Fitted line and 95% interval of the mean (normal approximation, n is in the thousands)."""
    n = len(x)
    residual_sd = np.sqrt(np.sum((y - np.polyval(fit, x)) ** 2) / (n - 2))
    sxx = np.sum((x - x.mean()) ** 2)
    half = 1.96 * residual_sd * np.sqrt(1 / n + (grid - x.mean()) ** 2 / sxx)
    line = np.polyval(fit, grid)
    return line, line - half, line + half


def name_similarity_vs_identity(result, stats, out_dir, name=NAME):
    """A: distribution of name similarity; B: distribution of sequence identity; C: density of the two against each
    other with both regressions and the 0.95 threshold. Writes the figure (svg, png) and source-data CSVs to out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    table = result["distance"]
    a = stats["analysis"]
    x = table["sequence_identity"].to_numpy()
    y_full = table["name_similarity"].to_numpy()
    y_free = result["distance_ignoring_alleles"]["name_similarity_without_alleles"].to_numpy()
    n = len(table)
    assert n == a["matches_analysed"]

    hist_a, hist_b = _histogram(y_full), _histogram(x)
    counts, x_edges, y_edges = np.histogram2d(x, y_full, bins=[BINS_X, BINS_Y], range=[[0, 1], [0, 1]])
    assert hist_a["matches"].sum() == hist_b["matches"].sum() == counts.sum() == n
    ix, iy = np.nonzero(counts)
    density = pd.DataFrame({"identity_low": x_edges[ix], "identity_high": x_edges[ix + 1],
                            "name_similarity_low": y_edges[iy], "name_similarity_high": y_edges[iy + 1],
                            "matches": counts[ix, iy].astype(int)})

    grid = np.linspace(x.min(), 1.0, 101)
    fit_rows = []
    for label, y, block in (("full_names", y_full, a["full_names"]), ("without_alleles", y_free, a["without_alleles"])):
        fit = [block["slope"], block["intercept"]]
        line, lower, upper = _fit_band(x, y, fit, grid)
        fit_rows.append(pd.DataFrame({"analysis": label, "sequence_identity": grid, "fit": line,
                                      "lower_95": lower, "upper_95": upper}))
    fits = pd.concat(fit_rows, ignore_index=True)

    fig = style.new_figure(WIDTH_MM, HEIGHT_MM)
    ax_a = style.axes_mm(fig, 20, 100, 65, 42)
    ax_b = style.axes_mm(fig, 112, 100, 64, 42)
    ax_c = style.axes_mm(fig, 20, 14, 142, 70)
    ax_bar = style.axes_mm(fig, 166, 14, 4, 70)

    for ax, hist, xlabel, letter in ((ax_a, hist_a, "Name similarity (normalised Levenshtein)", "A"),
                                     (ax_b, hist_b, "Sequence identity", "B")):
        ax.bar(hist["bin_low"], hist["matches"], width=hist["bin_high"] - hist["bin_low"], align="edge",
               color=style.GREY, linewidth=0)
        ax.set_xlim(0, 1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Best matches")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
        ax.text(-0.24, 1.05, letter, transform=ax.transAxes, fontsize=8, fontweight="bold")

    mesh = ax_c.pcolormesh(x_edges, y_edges, np.where(counts.T > 0, counts.T, np.nan), norm=LogNorm(vmin=1, vmax=counts.max()),
                           cmap="viridis", rasterized=True)
    colour_bar = fig.colorbar(mesh, cax=ax_bar)
    colour_bar.set_label("Best matches per bin (log scale)")
    for label, colour, linestyle, text in (("full_names", style.FULL, "-", "full names"),
                                           ("without_alleles", style.ALLELE_FREE, "--", "alleles ignored")):
        part = fits[fits["analysis"] == label]
        ax_c.fill_between(part["sequence_identity"], part["lower_95"], part["upper_95"], color=colour, alpha=0.35, linewidth=0)
        ax_c.plot(part["sequence_identity"], part["fit"], color=colour, linestyle=linestyle, linewidth=1.0,
                  label=f"{text}: y = {a[label]['slope']:.2f}x {'-' if a[label]['intercept'] < 0 else '+'} {abs(a[label]['intercept']):.2f}")
    ax_c.axvline(0.95, color="white", linestyle=":", linewidth=0.8)
    for label, colour, dy in (("full_names", style.FULL, -0.07), ("without_alleles", style.ALLELE_FREE, 0.05)):
        predicted = a[label]["predicted_similarity_at_threshold"]
        ax_c.plot([0.95], [predicted], marker="o", markersize=3.5, color=colour, markeredgecolor="white", markeredgewidth=0.4)
        ax_c.text(0.93, predicted + dy, f"{predicted:.2f}", ha="right", va="center", color="black", fontsize=7,
                  bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": colour, "linewidth": 0.6})
    ax_c.text(0.95, 1.015, "identity 0.95", color="black", ha="center", va="bottom", fontsize=7, clip_on=False)
    ax_c.text(0.60, 1.015, f"n = {n:,} best matches", color="black", ha="left", va="bottom", fontsize=7, clip_on=False)
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.set_xlabel("Sequence identity")
    ax_c.set_ylabel("Name similarity")
    ax_c.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, frameon=False, labelcolor="black", handlelength=1.6,
                columnspacing=1.0, borderaxespad=0.2)
    ax_c.text(-0.13, 1.1, "C", transform=ax_c.transAxes, fontsize=8, fontweight="bold")

    style.save_figure(fig, out_dir / name)
    for suffix, frame in (("A_histogram", hist_a), ("B_histogram", hist_b), ("C_density", density), ("C_fit", fits)):
        frame.to_csv(out_dir / f"source_data_{name}_{suffix}.csv", index=False)
    return fig
