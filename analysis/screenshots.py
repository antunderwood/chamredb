"""Terminal screenshots of `chamredb query`, rendered from the tool's own output (spec #34, ticket #42).

capture_query runs the package's query function with its printing redirected, so the text is exactly what the CLI prints.
render_screenshot writes an SVG made by rich (real emoji, real colours) and a PNG drawn by matplotlib in DejaVu Sans Mono with
the same colours; the emoji icons are replaced in the PNG by plain symbols of the same cell width because no bundled font has them.
"""
from io import StringIO
from pathlib import Path

from rich.cells import cell_len
from rich.console import Console
from rich.terminal_theme import MONOKAI
from rich.text import Text

from analysis import style

WIDTH = 100                 # terminal columns
FONT_PT = 7
CELL_MM = 0.60205 * FONT_PT * 25.4 / 72    # DejaVu Sans Mono advance
LINE_MM = 3.5
BACKGROUND = "#{:02x}{:02x}{:02x}".format(*MONOKAI.background_color)
FOREGROUND = "#{:02x}{:02x}{:02x}".format(*MONOKAI.foreground_color)
# emoji icons used by the CLI -> symbols DejaVu Sans Mono has, drawn in the same two cells
SUBSTITUTES = {"\U0001F9EC": "●", "\U0001F4C4": "▪", "\U0001F5C4": "▣", "\U0001F517": "∘", "➡": "→"}


def capture_query(identifier, database, graph_path=None, coverage=0.9, identity=0.9):
    """(the command, the markup strings the CLI prints) for `chamredb query -d DATABASE -i IDENTIFIER`."""
    import chamredb.query_graph as query_graph

    printed = []
    original = query_graph.rprint
    query_graph.rprint = printed.append
    try:
        query_graph.query_graph_single_id(identifier, database, coverage, identity, graph_path)
    finally:
        query_graph.rprint = original
    shown = f'"{identifier}"' if any(c in identifier for c in "() ") else identifier
    return f"chamredb query -d {database} -i {shown}", printed


def _lines(command, markups, emoji=True):
    console = Console(width=WIDTH, file=StringIO(), force_terminal=True, color_system="truecolor")
    out = [Text(f"$ {command}", style="bold")]
    for markup in markups:
        out.append(Text.from_markup(markup, emoji=emoji))
    segments = []
    for renderable in out:
        segments.extend(console.render_lines(renderable, console.options, pad=False))
    return segments


def plain_text(markups, width=WIDTH):
    console = Console(width=width, file=StringIO(), force_terminal=False)
    for markup in markups:
        console.print(markup)
    return console.file.getvalue()


def drawn_text(command, markups):
    """Every character the PNG draws (after the emoji substitution)."""
    return "".join(SUBSTITUTES.get(c, c) for line in _lines(command, markups) for segment in line for c in segment.text)


def _colour(style_):
    if style_ is not None and style_.color is not None:
        triplet = style_.color.get_truecolor(MONOKAI)
        return f"#{triplet.red:02x}{triplet.green:02x}{triplet.blue:02x}"
    return FOREGROUND


def render_screenshot(command, markups, out_dir, name):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = _lines(command, markups)
    height_mm = len(lines) * LINE_MM + 10
    width_mm = 165
    fig = style.new_figure(width_mm, height_mm)
    fig.patch.set_facecolor(BACKGROUND)
    ax = style.axes_mm(fig, 0, 0, width_mm, height_mm)
    ax.set_xlim(0, width_mm)
    ax.set_ylim(0, height_mm)
    ax.axis("off")
    for row, line in enumerate(lines):
        y = height_mm - 5 - row * LINE_MM
        x = 6.0
        for segment in line:
            weight = "bold" if segment.style is not None and segment.style.bold else "normal"
            colour = _colour(segment.style)
            run = ""
            for char in segment.text:
                if char in SUBSTITUTES or cell_len(char) != 1:
                    if run:
                        ax.text(x, y, run, family="DejaVu Sans Mono", fontsize=FONT_PT, color=colour, fontweight=weight, ha="left", va="top")
                        x += len(run) * CELL_MM
                        run = ""
                    ax.text(x, y, SUBSTITUTES.get(char, char), family="DejaVu Sans Mono", fontsize=FONT_PT, color=colour, fontweight=weight, ha="left", va="top")
                    x += cell_len(char) * CELL_MM
                else:
                    run += char
            if run:
                ax.text(x, y, run, family="DejaVu Sans Mono", fontsize=FONT_PT, color=colour, fontweight=weight, ha="left", va="top")
                x += len(run) * CELL_MM
    style.save_figure(fig, out_dir / name)

    console = Console(record=True, width=WIDTH, file=StringIO(), force_terminal=True, color_system="truecolor")
    console.print(Text(f"$ {command}", style="bold"))
    for markup in markups:
        console.print(markup)
    (out_dir / f"{name}.svg").write_text(console.export_svg(title=command, theme=MONOKAI))
    plain = Console(width=WIDTH, file=StringIO(), force_terminal=False)
    plain.print(f"$ {command}")
    for markup in markups:
        plain.print(markup)
    (out_dir / f"source_data_{name}.txt").write_text(plain.file.getvalue())
    return fig
