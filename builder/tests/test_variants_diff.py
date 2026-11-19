"""Release diff of the mutation tables."""
import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "builder" / "scripts"))
from functions import variants_diff  # noqa: E402


def table():
    rec = lambda rid, db: {"record_id": rid, "database": db, "native_notation": rid, "taxon_group": "Escherichia"}
    return {
        "records": [rec("a", "ncbi"), rec("b", "card")],
        "matches": [{"record_id": "a", "other_database": "card", "class": "same_mutation"},
                    {"record_id": "b", "other_database": "ncbi", "class": "same_mutation"}],
        "loci": {"Escherichia:gyrA:protein": {"anchor_database": "ncbi", "anchor_accession": "WP_1.1", "anchor_sha256": "x"}},
        "clusters": [{"members": ["a", "b"]}],
    }


def test_no_change():
    d = variants_diff.diff(table(), table())
    assert not d["added"] and not d["removed"] and not d["class_changes"] and not d["anchor_changes"]


def test_added_removed_class_and_anchor_changes():
    old, new = table(), copy.deepcopy(table())
    new["records"].append({"record_id": "c", "database": "pointfinder", "native_notation": "c", "taxon_group": "Escherichia"})
    new["records"] = [r for r in new["records"] if r["record_id"] != "b"]
    new["matches"][0]["class"] = "locus_not_covered"
    new["loci"]["Escherichia:gyrA:protein"]["anchor_accession"] = "WP_2.1"
    d = variants_diff.diff(old, new)
    assert d["added"] == ["c"] and d["removed"] == ["b"]
    assert d["class_changes"] == [(("a", "card"), "same_mutation", "locus_not_covered")]
    assert d["anchor_changes"] == [("Escherichia:gyrA:protein", "ncbi WP_1.1", "ncbi WP_2.1")]
    text = variants_diff.markdown(d)
    assert "Anchor changes (reported positions move at these loci) (1)" in text and "same_mutation -> locus_not_covered: 1" in text


def test_bundled_table_against_itself():
    data = json.loads((ROOT / "chamredb" / "data" / "variants.json").read_text())
    d = variants_diff.diff(data, data)
    assert d["records"][0] == d["records"][1] and not d["class_changes"]
