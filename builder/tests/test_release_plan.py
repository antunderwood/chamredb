"""Release planning (ticket #8, ADR 0005): tag, change detection and notes, on hand-made manifests (level 1)."""
import datetime
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "release"))
import plan_release as pr  # noqa: E402
import resolve_inputs as ri  # noqa: E402


def manifest(card="4.0.2", ncbi="AMRFinderPlus database 4.2, 2026-10-29.1", resfinder="eecf0aa"):
    f = lambda db, role, release, **kw: {"database": db, "role": role, "release": release, **kw}
    return {"files": [f("card", "sequences", card, card_release=card), f("card", "ontology", card, card_release=card),
                      f("ncbi", "catalogue", ncbi), f("ncbi", "proteins", ncbi),
                      f("resfinder", "phenotypes", resfinder), f("resfinder", "nucleotides", resfinder)]}


def inputs(card="4.0.2", ncbi="AMRFinderPlus database 4.2, 2026-10-29.1", resfinder="eecf0aa"):
    return {"card": {"release": card}, "ontology": {"release": card}, "ncbi": {"release": ncbi},
            "resfinder": {"release": f"resfinder_db {resfinder}", "commit": resfinder + "0" * 33}}


def test_tag_is_dated():
    assert pr.tag_for(datetime.date(2026, 11, 9)) == "graph-2026.11.09"


def test_releases_from_manifest_and_inputs_agree_when_nothing_changed():
    assert pr.releases_from_manifest(manifest()) == pr.releases_from_inputs(inputs())
    assert pr.releases_from_manifest(manifest()) == {"card": "4.0.2", "ncbi": "AMRFinderPlus database 4.2, 2026-10-29.1", "resfinder": "eecf0aa"}


def test_no_previous_release_means_build():
    build, reasons = pr.needs_build(None, inputs())
    assert build and reasons == ["no previous release"]


def test_unchanged_upstream_means_no_build():
    assert pr.needs_build(manifest(), inputs()) == (False, [])


@pytest.mark.parametrize("change,expected", [
    ({"card": "4.0.3"}, ["card: 4.0.2 -> 4.0.3"]),
    ({"resfinder": "abcdef1"}, ["resfinder: eecf0aa -> abcdef1"]),
    ({"ncbi": "AMRFinderPlus database 4.2, 2026-11-30.1"}, ["ncbi: AMRFinderPlus database 4.2, 2026-10-29.1 -> AMRFinderPlus database 4.2, 2026-11-30.1"]),
])
def test_any_changed_database_triggers_a_build_and_is_named(change, expected):
    build, reasons = pr.needs_build(manifest(), inputs(**change))
    assert build and reasons == expected


def test_release_notes_list_releases_and_embed_the_report():
    notes = pr.release_notes("graph-2026.11.09", manifest(), "# Graph diff\n- Nodes gained: 5\n")
    for line in ("graph-2026.11.09", "CARD 4.0.2", "AMRFinderPlus database 4.2, 2026-10-29.1", "eecf0aa", "Nodes gained: 5"):
        assert line in notes


def test_card_version_is_read_from_card_json():
    assert ri.parse_card_version(json.dumps({"_version": "4.0.2", "1": {}}).encode()) == "4.0.2"
    with pytest.raises(ValueError, match="_version"):
        ri.parse_card_version(b"{}")


def test_ncbi_directory_is_the_newest_major_directory_that_holds_the_version():
    listing = '<a href="3.9/">3.9/</a><a href="4.0/">4.0/</a><a href="4.2/">4.2/</a><a href="changelogs/">c</a>'
    has = {("4.0", "2026-10-29.1"): False, ("4.2", "2026-10-29.1"): True}
    assert ri.pick_ncbi_dir("2026-10-29.1", listing, lambda d, v: has.get((d, v), False)) == "4.2"
    with pytest.raises(ValueError, match="no database directory"):
        ri.pick_ncbi_dir("2099-01-01.1", listing, lambda d, v: False)
