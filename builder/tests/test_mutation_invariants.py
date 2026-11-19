"""The bundled variants.json passes the checks, and the tables in builder/variants are consistent."""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "builder" / "scripts"))
from functions import mutation_compare as mc  # noqa: E402
from functions import mutation_invariants  # noqa: E402

BUNDLED = ROOT / "chamredb" / "data" / "variants.json"
TABLES = ROOT / "builder" / "variants"


@pytest.fixture(scope="module")
def data():
    return json.loads(BUNDLED.read_text())


def test_bundled_variants_hold_the_invariants(data):
    assert mutation_invariants.check_variants(data) == {}


def test_every_record_matched_against_both_other_databases(data):
    assert len(data["matches"]) == 2 * len(data["records"])


def test_taxa_table_is_unique_and_complete():
    rows = mc.read_table(TABLES / "taxa.tsv")
    keys = [(r["database"], r["label"]) for r in rows]
    assert len(keys) == len(set(keys))
    families = {}
    for r in rows:   # a group has one family and one source wherever it appears
        families.setdefault(r["group"], set()).add((r["group_family"], r["group_source"]))
    assert all(len(v) == 1 for v in families.values()), {k: v for k, v in families.items() if len(v) > 1}
    # NCBI groups are its own AMRFinderPlus organisms; PointFinder groups are AMRFinderPlus organisms or its species
    amrfinder = {r["group"] for r in rows if r["database"] == "ncbi"}
    assert all(r["group"] == r["label"] for r in rows if r["database"] == "ncbi")
    assert all(r["group"] in amrfinder or r["group"] == r["label"] for r in rows if r["database"] == "pointfinder")


def test_drug_class_table_is_unique():
    rows = mc.read_table(TABLES / "drug_classes.tsv")
    keys = [(r["database"], r["term"]) for r in rows]
    assert len(keys) == len(set(keys))


def test_reviewed_equivalences_name_existing_records(data):
    ids = {r["record_id"] for r in data["records"]}
    for row in mc.read_table(TABLES / "reviewed_equivalences.tsv"):
        assert row["record_a"] in ids and row["record_b"] in ids, row
        assert row["decision"] in ("same", "different"), row
