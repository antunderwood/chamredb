"""Figures of the README analysis (spec #34). Display transforms only (binning, a fitted line and its 95% interval);
every statistic comes from the result of analysis.matches.assess and analysis.run.statistics."""
import textwrap
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


OVERLAP_NAME = "figure_overlap_upset"
DATABASES = ["card", "ncbi", "resfinder"]
DATABASE_LABELS = {"card": "CARD", "ncbi": "NCBI", "resfinder": "ResFinder"}


def overlap_upset(group_stats, out_dir, name=OVERLAP_NAME):
    """UpSet figure of the equivalence groups (proteins joined by reciprocal best hits) by the databases they contain.
    group_stats is the `groups` block of statistics.json. Writes the figure (svg, png) and two source-data CSVs."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for profile, groups in group_stats["by_databases"].items():
        members = profile.split("+")
        rows.append({"databases": profile, **{d: int(d in members) for d in DATABASES},
                     "groups": groups, "proteins": group_stats["proteins_by_databases"][profile]})
    inter = pd.DataFrame(rows).sort_values(["groups", "databases"], ascending=[False, True]).reset_index(drop=True)
    sets = pd.DataFrame({"database": DATABASES,
                         "groups": [int(inter.loc[inter[d] == 1, "groups"].sum()) for d in DATABASES]})
    total, proteins = group_stats["total"], int(inter["proteins"].sum())
    assert inter["groups"].sum() == total

    fig = style.new_figure(180, 105)
    ax_bars = style.axes_mm(fig, 62, 50, 114, 48)
    ax_matrix = style.axes_mm(fig, 62, 12, 114, 30)
    ax_sets = style.axes_mm(fig, 10, 12, 32, 30)
    positions = np.arange(len(inter))

    ax_bars.bar(positions, inter["groups"], width=0.6, color=style.GREY, linewidth=0)
    for x, value in zip(positions, inter["groups"]):
        ax_bars.text(x, value, f"{value:,}", ha="center", va="bottom", fontsize=7)
    ax_bars.set_xlim(-0.6, len(inter) - 0.4)
    ax_bars.set_ylim(0, inter["groups"].max() * 1.12)
    ax_bars.set_xticks([])
    ax_bars.spines["bottom"].set_visible(False)
    ax_bars.set_ylabel("Equivalence groups")
    ax_bars.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    ax_bars.text(0.99, 0.97, f"{total:,} groups from {proteins:,} proteins", transform=ax_bars.transAxes,
                 ha="right", va="top", fontsize=7)

    ax_matrix.set_xlim(ax_bars.get_xlim())
    ax_matrix.set_ylim(-0.6, len(DATABASES) - 0.4)
    ax_matrix.invert_yaxis()
    for x, row in inter.iterrows():
        members = [i for i, d in enumerate(DATABASES) if row[d]]
        ax_matrix.plot([x, x], [min(members), max(members)], color="black", linewidth=1.0, zorder=1)
        for i, d in enumerate(DATABASES):
            ax_matrix.scatter([x], [i], s=22, color="black" if row[d] else "#DDDDDD", zorder=2, linewidths=0)
    ax_matrix.set_yticks(range(len(DATABASES)))
    ax_matrix.set_yticklabels([DATABASE_LABELS[d] for d in DATABASES])
    ax_matrix.set_xticks([])
    for spine in ("left", "bottom"):
        ax_matrix.spines[spine].set_visible(False)
    ax_matrix.tick_params(axis="y", length=0)

    ax_sets.barh(range(len(DATABASES)), sets["groups"], color=[style.DATABASE_COLOURS[d] for d in DATABASES], height=0.6, linewidth=0)
    for i, value in enumerate(sets["groups"]):
        ax_sets.annotate(f"{value:,}", xy=(value, i), xytext=(-3, 0), textcoords="offset points", ha="right", va="center", fontsize=7)
    ax_sets.set_xlim(sets["groups"].max() * 1.45, 0)
    ax_sets.set_ylim(-0.6, len(DATABASES) - 0.4)
    ax_sets.invert_yaxis()
    ax_sets.set_yticks([])
    ax_sets.spines["left"].set_visible(False)
    ax_sets.spines["right"].set_visible(True)
    ax_sets.yaxis.tick_right()
    ax_sets.set_xlabel("Groups containing\nthe database")
    ax_sets.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))

    fig.text(8 / 180, 98 / 105, "\n".join(textwrap.wrap(
        "Equivalence groups are proteins joined by reciprocal best hits across the three databases. "
        "Each bar counts the groups that contain exactly the databases marked by the dots below it.", 36)),
        ha="left", va="top", fontsize=7)
    style.save_figure(fig, out_dir / name)
    inter.to_csv(out_dir / f"source_data_{name}_intersections.csv", index=False)
    sets.to_csv(out_dir / f"source_data_{name}_set_sizes.csv", index=False)
    return fig
