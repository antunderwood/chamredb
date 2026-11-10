# Analysis harness (tickets #35, #37, #38)

`make test-analysis` runs `test_matches.py`, `test_cleaning.py` and `test_run.py` (about 3 s; the last also runs the command twice on the bundled release graph). It needs the optional extra: `pip install -e ".[analysis]" pytest`.

## Oracles

1. **Hand graph (level 1).** 12 nodes, 13 best matches, one skipped `blaPAO` match. Name similarities are exact fractions worked
   out by hand in the test (qnrb2 vs qnrb20 = 5/6, vana vs vanhax = 2/3, the rest identical); the regression is checked against
   closed-form least squares with `fractions`, not numpy. One case fixes the quirk below.
2. **Frozen output of the unmodified script (level 6).** `scripts/assess_matches.py` run unchanged on the committed 176-node
   subset of the 2021 graph (`chamredb/tests/test_data/graph_2021_subset.json`): the two CSVs in `expected/` and the statistics
   in the test. It proves "unchanged", not "correct". Environment of the freeze (2026-11-09): Python 3.11.17, networkx 3.4.2,
   pandas 3.1.1, numpy 2.4.6, scipy 1.17.1, strsimpy 0.4.1, seaborn 0.13.2, matplotlib 3.11.2.
   The script crashes on networkx 3.6 (`node_link_graph` default key `edges`, as found in #5), so 3.4.2 was used.

## Mutation check (the harness must go red)

Run 2026-11-09 in the working tree, then restored:

| Mutation in `analysis/matches.py` | Tests that fail |
| --- | --- |
| name cleaning keeps the `bla` prefix | hand graph (best matches, regression) and the three frozen-script tests |
| best match ignores identity | the three frozen-script tests |
| no `blaPAO` skip rule | the two hand-graph tests |

Added in ticket #37 (run 2026-11-09): pattern two group back to pattern one fails the `dfrA6 from Proteus mirabilis` case and the
bundled-graph invariants; comparing a name with itself instead of its match fails the differential and invariant tests; allowing
digits at the start of the allele pattern fails the frozen-script tests. `1 - distance` instead of `similarity` survives because
strsimpy's normalised distance is already distance over longest length (an equivalent mutant, not a gap).

Added in ticket #38 (run 2026-11-09): `bla` removed everywhere again fails `test_bla_is_removed_only_as_a_prefix`; the allele-free
threshold from the full-name fit fails `test_allele_free_threshold_comes_from_its_own_regression`; RBH pairs counted per direction fails
`test_statistics_on_the_hand_graph`.

## Decisions of ADR 0006 and `legacy=True`

The frozen-script tests (level 6) call the library with `legacy=True`, which keeps the two rules ADR 0006 changed: `bla` removed
anywhere, and the allele-free threshold from the full-name regression. Everything else is the same code for both, so the frozen
CSVs still prove that the library equals the script apart from those two named rules. The default (`legacy=False`) applies ADR 0006.

## Behaviour kept on purpose by `legacy=True` (defects, now changed by default except where noted)

- **Best match depends on edge order.** A later target replaces the best only if its coverage *and* its identity both exceed the
  current best (`best_target`), so a higher-identity target with lower coverage is ignored. Hand-graph case: c4 keeps vanHAX
  (identity 0.9, coverage 0.95) over blaTEM-1 (0.95, 0.9).
- **The allele-free threshold uses the wrong regression** (changed by ADR 0006). `predicted_similarity_without_alleles` is computed from the full-name
  fit, so both thresholds are the same number (test asserts equality). This is probably the source of the README's 0.69 and 0.86.
- ~~Second CARD species pattern~~ fixed in ticket #37: it took its group from the first pattern (would raise `AttributeError`)
  and its epithet was one letter, so it never matched a real name. In the 2026-11-09 graph it now applies to one CARD name
  (`dfrA6 from Proteus mirabilis` becomes `dfra6`); patterns one and three apply to 3 and 46 names.
- **`bla` is removed wherever it occurs, not only as a prefix** (README says prefix): 6 names are affected, for example
  `cblA` becomes `c`. Changed by ADR 0006 (prefix only).
- **Skip rule:** the code skips only nodes named exactly `blaPAO`; the README also names blaPDC (decision #36).
- A node whose targets all have coverage 0 would crash the script; the library skips it (cannot occur on a valid graph: coverage > 0).

## Facts recorded, not tested

Unmodified script on the full 2021 graph (11,718 nodes, 33,632 edges, 3.6 s): 22,276 best matches, 955 `blaPAO` skipped,
name similarity = 1.3431 x identity - 0.3989, predicted similarity at identity 0.95 = **0.877** (the README says 0.69),
19,618 matches with identity >= 0.95, 2,360 rows in the low-name-similarity table and 369 without alleles (the committed
CSVs in `scripts/` have 2,715 and 389). The committed README numbers, CSVs and figures therefore do not come from the shipped 2021
graph either.

## Figure style (ticket #39)

`analysis/style.py` is the shared style: sizes in millimetres (`new_figure(w_mm, h_mm)`, `axes_mm`), fonts of at least 7 pt, the
Okabe-Ito palette (one colour per database), editable text (`svg.fonttype none`, `pdf.fonttype 42`), and `save_figure`, which writes SVG
and a 300 dpi PNG with no date or software metadata so re-runs are byte-identical. `analysis/figures.py` draws the figures and writes a
`source_data_*` CSV for every panel beside them (`docs/images/`). Tests check sizes, fonts, palette, that the source data sum to the
statistics, that numbers written on the figure equal the statistics, and byte-identity across renders; the picture itself is reviewed
by the maintainer. Mutation (run 2026-11-09): a fit line with a 10% steeper slope fails `test_source_data_equals_the_statistics`.

Byte-identical outputs are guaranteed between runs in one environment (tested). Another platform's numpy gives regression coefficients that
agree to about 1e-13 but not bit for bit (found when CI disagreed with a Mac run), so tests that compare with the committed
`docs/statistics.json` use a tolerance of 1e-9 for floats and exact equality for counts.
