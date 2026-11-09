"""CARD 4.0.2 aro.obo has `property_value: <relation>` lines with no value (1,402 of them), which pronto rejects.
The builder drops exactly those lines before parsing (ticket #12); the downloaded file stays untouched.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "functions"))
import obo_cleaning as oc  # noqa: E402

OBO = """[Term]
id: ARO:1
name: a
property_value: classified_as_drug_class
is_a: ARO:2 ! b
property_value: has_part ARO:3
property_value: note "text" xsd:string
property_value: classified_as_antibiotic
"""


def test_drops_only_valueless_property_value_lines_and_counts_them():
    cleaned, n = oc.strip_valueless_property_values(OBO)
    assert n == 2
    assert "property_value: classified_as_drug_class" not in cleaned
    assert "property_value: classified_as_antibiotic" not in cleaned
    for kept in ("id: ARO:1", "is_a: ARO:2 ! b", "property_value: has_part ARO:3", 'property_value: note "text" xsd:string'):
        assert kept in cleaned
    assert cleaned.count("\n") == OBO.count("\n") - 2


def test_clean_input_is_unchanged():
    text = "[Term]\nid: ARO:1\nname: a\n"
    assert oc.strip_valueless_property_values(text) == (text, 0)


def test_refuses_an_implausible_number_of_removals():
    many = "[Term]\n" + "property_value: x\n" * 50
    with pytest.raises(ValueError, match="more than"):
        oc.strip_valueless_property_values(many, max_removed=10)
