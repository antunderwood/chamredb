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


HERO_NAME = "figure_real_families"
_ORDER = ["card", "ncbi", "resfinder"]


def real_families(graph, out_dir, name=HERO_NAME):
    """Six real families chosen by the rules of analysis.families, drawn with the databases in three columns: node colour
    is the database, solid lines are reciprocal best hits, dashed arrows are one-way hits, identity is written on the edges
    (when a panel has more than four edges, one note states that every edge has identity 1.00). Writes the figure and two
    source-data CSVs (every drawn node and edge)."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch
    from analysis import families

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chosen = families.select_families(graph)
    assert len(chosen) == 6

    fig = style.new_figure(180, 100)
    lefts, bottoms = [6, 66, 126], [56, 12]
    node_rows, edge_rows = [], []
    for index, family in enumerate(chosen):
        ax = style.axes_mm(fig, lefts[index % 3], bottoms[index // 3], 52, 30)
        ax.axis("off")
        by_db = {d: [n for n in family["nodes"] if n.split(":")[0] == d] for d in _ORDER}
        tallest = max(len(v) for v in by_db.values())
        ax.set_xlim(-0.55, 3.0)
        ax.set_ylim(-2.5, 2.5)   # the same scale in every panel
        position = {}
        for column, database in enumerate(_ORDER):
            members = by_db[database]
            for row, node in enumerate(members):
                position[node] = (column, (len(members) - 1) / 2 - row)
            ax.text(column, 2.45, style_label(database), ha="center", va="top", fontsize=7,
                    color=style.DATABASE_COLOURS[database], fontweight="bold")
        letter = "ABCDEF"[index]
        top = ax.get_xaxis_transform()   # x in data units, y as a fraction of the axes
        ax.text(-0.55, 1.02, letter, fontsize=8, fontweight="bold", ha="left", va="bottom", clip_on=False, transform=top)
        ax.text(-0.3, 1.02, family["title"], fontsize=7, ha="left", va="bottom", clip_on=False, transform=top)

        label_edges = len(family["edges"]) <= 4
        for source, target in family["edges"]:
            data = graph.get_edge_data(source, target)
            (x0, y0), (x1, y1) = position[source], position[target]
            if abs(x1 - x0) < 2:
                rad = 0.0
            elif y1 > y0:
                rad = -0.25      # bow upwards for targets above the source
            elif y1 < y0:
                rad = 0.25       # bow downwards for targets below the source
            else:
                rad = 0.9        # same height: a wide arc below the middle column
            owh = data["type"] == "OWH"
            ax.add_patch(FancyArrowPatch(
                (x0, y0), (x1, y1), arrowstyle="-|>" if owh else "-", mutation_scale=6, shrinkA=4, shrinkB=4,
                linestyle="--" if owh else "-", color=style.GREY if owh else "black", linewidth=0.8,
                connectionstyle=f"arc3,rad={rad}", zorder=1))
            if label_edges:
                mx, my = (x0 + x1) / 2, (y0 + y1) / 2
                dx, dy = x1 - x0, y1 - y0
                mx, my = mx + 0.5 * rad * dy, my - 0.5 * rad * dx   # apex of the arc
                ax.text(mx, my - (0.12 if rad == 0 else 0.1), f"{data['identity']:.2f}", ha="center", va="top", fontsize=7,
                        bbox={"boxstyle": "round,pad=0.1", "facecolor": "white", "edgecolor": "none", "alpha": 0.9}, zorder=3)
            edge_rows.append({"panel": letter, "key": family["key"], "source": source, "target": target,
                              "type": data["type"], "identity": data["identity"], "coverage": data["coverage"]})
        if not label_edges:
            ax.text(1.0, -2.45, f"identity 1.00 on all {len(family['edges'])} edges", ha="center", va="bottom", fontsize=7)
        for node, (x, y) in position.items():
            database = node.split(":")[0]
            ax.scatter([x], [y], s=60, color=style.DATABASE_COLOURS[database], edgecolors="white", linewidths=0.6, zorder=4)
            stacked = database == "resfinder" and len(by_db["resfinder"]) > 1   # names to the right of stacked nodes
            ax.text(x + (0.12 if stacked else 0), y + (0 if stacked else 0.2), graph.nodes[node]["name"],
                    ha="left" if stacked else "center", va="center" if stacked else "bottom", fontsize=7, zorder=5,
                    bbox={"boxstyle": "round,pad=0.05", "facecolor": "white", "edgecolor": "none", "alpha": 0.85})
            node_rows.append({"panel": letter, "key": family["key"], "database": database, "node_id": node,
                              "name": graph.nodes[node]["name"]})

    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6, color=style.DATABASE_COLOURS[d], label=style_label(d)) for d in _ORDER]
    handles += [Line2D([0], [0], color="black", linewidth=0.8, label="reciprocal best hit"),
                Line2D([0], [0], color=style.GREY, linewidth=0.8, linestyle="--", label="one-way hit (arrow: from, to)")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.0), handlelength=1.8, columnspacing=1.2)

    style.save_figure(fig, out_dir / name)
    pd.DataFrame(node_rows).to_csv(out_dir / f"source_data_{name}_nodes.csv", index=False)
    pd.DataFrame(edge_rows).to_csv(out_dir / f"source_data_{name}_edges.csv", index=False)
    return fig


def style_label(database):
    return DATABASE_LABELS[database]


DIST_NAME = "figure_name_distributions"
GALLERY_NAME = "figure_name_disagreements"


def name_distributions(result, stats, out_dir, name=DIST_NAME):
    """Violin plots of the matches with identity at or above the threshold: A and B the name similarity (full names; alleles
    ignored), C and D the absolute difference from each analysis's own regression prediction. Writes the figure, a summary CSV
    (n, minimum, 2.5th and 97.5th percentiles, quartiles, median, maximum per panel) and a histogram CSV."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    a = stats["analysis"]
    high, high_free = result["high_identity"], result["high_identity_without_alleles"]
    panels = [
        ("A", "full_names", "Name similarity, full names", high["name_similarity"].to_numpy(), "Name similarity"),
        ("B", "without_alleles", "Name similarity, alleles ignored", high_free["name_similarity_without_alleles"].to_numpy(), "Name similarity"),
        ("C", "full_names", "Difference from the prediction, full names", high["similarity_difference"].to_numpy(), "|observed - predicted|"),
        ("D", "without_alleles", "Difference from the prediction, alleles ignored", high_free["similarity_difference"].to_numpy(), "|observed - predicted|"),
    ]
    fig = style.new_figure(180, 115)
    lefts, bottoms = [22, 108], [66, 14]
    summary_rows, hist_rows = [], []
    for index, (letter, block, title, values, ylabel) in enumerate(panels):
        ax = style.axes_mm(fig, lefts[index % 2], bottoms[index // 2], 62, 36)
        n = len(values)
        assert n == a[block]["matches_at_or_above_threshold"]
        parts = ax.violinplot(values, positions=[0], widths=0.9, showextrema=False, showmedians=False)
        for body in parts["bodies"]:
            body.set_facecolor(style.GREY)
            body.set_alpha(0.55)
            body.set_linewidth(0)
        percentiles = np.percentile(values, [2.5, 25, 50, 75, 97.5])
        ax.plot([0, 0], [percentiles[0], percentiles[4]], color="black", linewidth=1.0)
        ax.plot([-0.12, 0.12], [percentiles[2]] * 2, color=style.FULL, linewidth=1.6)
        ax.set_xlim(-0.6, 0.6)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
        ax.set_ylabel(ylabel)
        top = ax.get_xaxis_transform()
        ax.text(-0.6, 1.06, letter, fontsize=8, fontweight="bold", ha="left", va="bottom", clip_on=False, transform=top)
        ax.text(-0.5, 1.06, title, fontsize=7, ha="left", va="bottom", clip_on=False, transform=top)
        notes = [f"n = {n:,}", f"median {percentiles[2]:.2f}"]
        if letter in "AB":
            predicted = a[block]["predicted_similarity_at_threshold"]
            notes.append(f"{a[block]['low_name_similarity_rows'] / n:.1%} below {predicted:.2f}")
        else:
            notes.append(f"97.5% within {percentiles[4]:.2f}")
        ax.text(0.58, 0.97, "\n".join(notes), ha="right", va="top", fontsize=7, clip_on=False)
        summary_rows.append({"panel": letter, "analysis": block, "quantity": ylabel, "n": n, "min": values.min(),
                             "p2_5": percentiles[0], "p25": percentiles[1], "median": percentiles[2], "p75": percentiles[3],
                             "p97_5": percentiles[4], "max": values.max()})
        counts, edges = np.histogram(values, bins=20, range=(0, 1))
        assert counts.sum() == n
        hist_rows.append(pd.DataFrame({"panel": letter, "bin_low": edges[:-1], "bin_high": edges[1:], "matches": counts}))
    fig.text(6 / 180, 2 / 115, f"Matches with sequence identity of at least {stats['analysis']['identity_threshold']}; "
             "violin: kernel density over the observed range, black line: 2.5th to 97.5th percentile, orange: median.",
             ha="left", va="bottom", fontsize=7)
    style.save_figure(fig, out_dir / name)
    pd.DataFrame(summary_rows).to_csv(out_dir / f"source_data_{name}_summary.csv", index=False)
    pd.concat(hist_rows, ignore_index=True).to_csv(out_dir / f"source_data_{name}_histogram.csv", index=False)
    return fig


def name_disagreement_gallery(table, out_dir, name=GALLERY_NAME, per_pair=4):
    """A table figure of the most striking name disagreements (analysis.families.select_disagreements): identical or
    near-identical proteins whose names have nothing in common once alleles are ignored. Illustrative, not representative."""
    from analysis import families

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chosen = families.select_disagreements(table, per_pair=per_pair)
    everything = families.select_disagreements(table, per_pair=10 ** 9)
    eligible = len(everything)
    tied = int((everything["name_similarity_without_alleles"] == 0).sum())

    fig = style.new_figure(180, 95)
    width, height = fig.get_size_inches() * 25.4
    columns = [8, 58, 100, 140, 154]       # x in mm: databases, name in the first, name in the second, identity, similarity
    y = height - 6
    for x, header in zip(columns, ["Databases", "Name in the first", "Name in the second", "Identity", "Name similarity"]):
        fig.text(x / width, y / height, header, fontsize=7, fontweight="bold", ha="left", va="top")
    fig.add_artist(plt_line(8 / width, (y - 5) / height, 176 / width))
    y -= 7
    previous = None
    for r in chosen.itertuples():
        if r.database_pair != previous:
            previous = r.database_pair
        fig.text(columns[0] / width, y / height, f"{DATABASE_LABELS[r.database_a]} and {DATABASE_LABELS[r.database_b]}", fontsize=7, ha="left", va="top")
        fig.text(columns[1] / width, y / height, r.name_a, fontsize=7, ha="left", va="top", color=style.DATABASE_COLOURS[r.database_a], fontweight="bold")
        fig.text(columns[2] / width, y / height, r.name_b, fontsize=7, ha="left", va="top", color=style.DATABASE_COLOURS[r.database_b], fontweight="bold")
        fig.text(columns[3] / width, y / height, f"{r.sequence_identity:.2f}", fontsize=7, ha="left", va="top")
        fig.text(columns[4] / width, y / height, f"{r.name_similarity_without_alleles:.2f}", fontsize=7, ha="left", va="top")
        y -= 4.9
    fig.text(8 / width, 8 / height,
             f"First {per_pair} for each pair of databases among {eligible:,} best matches with identity of at least {families.MIN_IDENTICAL} and names of at most "
             f"{families.GALLERY_NAME_LENGTH} characters,\nordered by name similarity (alleles ignored, lowest first) then by name; {tied:,} of them have similarity 0, "
             "so these examples are illustrative, not representative.\nEvery low name similarity at identity of at least 0.95: docs/analysis/low_name_similarity_without_alleles.csv.",
             fontsize=7, ha="left", va="bottom")
    style.save_figure(fig, out_dir / name)
    chosen.to_csv(out_dir / f"source_data_{name}_rows.csv", index=False)
    pd.DataFrame([{"eligible_matches": eligible, "with_similarity_0": tied}]).to_csv(out_dir / f"source_data_{name}_eligible.csv", index=False)
    return fig


def plt_line(x0, y, x1):
    from matplotlib.lines import Line2D
    return Line2D([x0, x1], [y, y], color="black", linewidth=0.5)
