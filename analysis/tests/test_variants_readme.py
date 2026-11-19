"""The mutation numbers in the README match the bundled variants.json (headline numbers are recounted here)."""
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pytest

from analysis import readme, variants

ROOT = Path(__file__).resolve().parents[2]
README = (ROOT / "README.md").read_text()
DATA = json.loads((ROOT / "chamredb" / "data" / "variants.json").read_text())
VSTATS = json.loads((ROOT / "docs" / "variants" / "statistics.json").read_text())


def markers():
    return re.findall(r"<!-- vstat:([\w.]+) -->(.*?)<!-- /vstat -->", README)


def test_statistics_file_is_current():
    assert json.loads(json.dumps(variants.statistics(DATA), sort_keys=True)) == VSTATS


def test_every_marker_is_rendered_from_the_statistics():
    found = markers()
    assert found, "no vstat markers in the README"
    for key, value in found:
        assert value == readme.format_vstat(VSTATS, key), key


def independent():
    records = DATA["records"]
    dbs = Counter(r["database"] for r in records)
    groups = defaultdict(set)
    for r in records:
        groups[r["taxon_group"]].add(r["database"])
    card = [r for r in records if r["database"] == "card"]
    checks = Counter(c["ref_check"] for r in card for c in r["components"])
    return {
        "records": f"{len(records):,}", "ncbi": f"{dbs['ncbi']:,}", "card": f"{dbs['card']:,}", "pointfinder": f"{dbs['pointfinder']:,}",
        "combinations": f"{sum(1 for r in records if len(r['components']) > 1):,}",
        "groups": f"{len(groups):,}", "groups_two": f"{sum(1 for v in groups.values() if len(v) > 1):,}",
        "groups_three": f"{sum(1 for v in groups.values() if len(v) == 3):,}",
        "card_tb_share": f"{sum(1 for r in card if r['taxon_group'] == 'mycobacterium_tuberculosis') / len(card):.0%}",
        "self_card": f"{checks['mismatch'] / (checks['ok'] + checks['mismatch']):.1%}", "self_card_n": f"{checks['mismatch']:,}",
        "release": DATA["data_release"],
    }


@pytest.mark.parametrize("key", sorted(independent()))
def test_headline_numbers_against_an_independent_count(key):
    values = {k: v for k, v in markers()}
    assert values[key] == independent()[key]
