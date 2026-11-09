"""Best matches, name similarity and regression for the graph assessment (spec #34).

A port of the logic of scripts/assess_matches.py as functions on a graph read from a path, with no import of package
internals. Ticket #35 was a characterisation; `legacy=True` still reproduces the unmodified script. The two rules changed
by ADR 0006 (bla removed only as a prefix; the allele-free analysis uses its own regression) apply unless legacy=True.
Defects once listed in
analysis/tests/README.md (fixed later: the 'bla' skip rule, the regression used for the allele-free threshold; the second
CARD species pattern was fixed in ticket #37).
"""
import json
import re

import numpy as np
import pandas as pd
from networkx.readwrite import json_graph
from strsimpy import NormalizedLevenshtein

IDENTITY_THRESHOLD = 0.95
_levenshtein = NormalizedLevenshtein()


def read_graph(path):
    with open(path) as fh:
        return json_graph.node_link_graph(json.load(fh), edges="links")


def target_info(graph, node_id):
    """{database: {id: {'node': attributes, 'edge_data': attributes}}} for the out-edges of a node."""
    targets = {}
    for edge in graph.edges(node_id):
        database = edge[1].split(":")[0]
        identifier = "".join(edge[1].split(":")[1:])
        targets.setdefault(database, {})[identifier] = {
            "node": graph.nodes[edge[1]], "edge_data": graph.get_edge_data(edge[0], edge[1]), "node_id": edge[1]}
    return targets


def best_target(targets):
    """The existing rule with both thresholds 0: a later target replaces the best only if its coverage AND its
    identity both exceed the current best, so the result depends on edge order."""
    best, best_identity, best_coverage = None, 0, 0
    for target in targets.values():
        coverage, identity = target["edge_data"]["coverage"], target["edge_data"]["identity"]
        if coverage >= 0 and coverage > best_coverage:
            if identity >= 0 and identity > best_identity:
                best, best_coverage, best_identity = target, coverage, identity
    return best


def clean_name(name, database, remove_alleles=False, legacy=False):
    species_pattern_one = r"^[A-Z][a-z]+ [a-z]+ (.+) conferring resistance to .+$"
    species_pattern_two = r"^(.+) from [A-Z][a-z]+ [a-z]+$"  # ticket #37: the epithet was [a-z] (one letter), so it never matched
    species_pattern_three = r"^[A-Z][a-z]+ [a-z]+ (.+)$"
    if database == "card" and re.search(species_pattern_one, name):
        cleaned = re.search(species_pattern_one, name).group(1)
    elif database == "card" and re.search(species_pattern_two, name):
        cleaned = re.search(species_pattern_two, name).group(1)  # ticket #37: used pattern one
    elif database == "card" and re.search(species_pattern_three, name):
        cleaned = re.search(species_pattern_three, name).group(1)
    else:
        cleaned = name
    cleaned = cleaned.lower()
    if legacy:
        cleaned = cleaned.replace("bla", "")  # the script removed bla anywhere in the name
    elif cleaned.startswith("bla"):
        cleaned = cleaned[3:]                 # ADR 0006: only a leading prefix
    cleaned = cleaned.replace("(", "").replace(")", "").replace("-", "")
    if remove_alleles:
        allele_pattern = r"^([A-Za-z]+)[0-9]+[A-Za-z0-9.]*$"
        if re.search(allele_pattern, cleaned):
            cleaned = re.search(allele_pattern, cleaned).group(1)
    return cleaned


def best_matches(graph, legacy=False):
    """(table, number skipped): one row per Node and other database that has an edge, with the best match."""
    rows, skipped = [], 0
    for source_node in graph.nodes():
        targets = target_info(graph, source_node)
        for database in targets:
            match = best_target(targets[database])
            if match is None:
                continue
            source_name = graph.nodes[source_node]["name"]
            source_database = graph.nodes[source_node]["database"]
            target_name = match["node"]["name"]
            target_database = match["node"]["database"]
            if source_name == "blaPAO" or target_name == "blaPAO":  # defect kept: README also skips blaPDC (ticket #36)
                skipped += 1
                continue
            cleaned_source = clean_name(source_name, source_database, legacy=legacy)
            cleaned_target = clean_name(target_name, target_database, legacy=legacy)
            rows.append({
                "source_name": source_name, "cleaned_source_name": cleaned_source, "source_database": source_database,
                "target_name": target_name, "cleaned_target_name": cleaned_target, "target_database": target_database,
                "sequence_identity": match["edge_data"]["identity"],
                "name_similarity": _levenshtein.similarity(cleaned_source, cleaned_target),
                "source_id": source_node, "target_id": match["node_id"],
            })
    return pd.DataFrame(rows), skipped


def _fitted(table, column, fit):
    return table.assign(
        predicted_similarity=lambda x: np.polyval(fit, x["sequence_identity"]),
        similarity_difference=lambda x: abs(x[column] - x["predicted_similarity"]))


def assess(graph, legacy=False):
    """Everything the script computed: tables, both regressions and the two tables of low name similarity."""
    distance, skipped = best_matches(graph, legacy=legacy)
    fit = np.polyfit(distance["sequence_identity"], distance["name_similarity"], 1)
    fitted = _fitted(distance, "name_similarity", fit)
    predicted = np.polyval(fit, IDENTITY_THRESHOLD)
    high = fitted.query("sequence_identity >= @IDENTITY_THRESHOLD")
    low = fitted.query("sequence_identity >= @IDENTITY_THRESHOLD and name_similarity <= @predicted")

    ignoring = distance.copy()
    for kind in ("source", "target"):
        ignoring[f"cleaner_{kind}_name"] = ignoring.apply(
            lambda x: clean_name(x[f"{kind}_name"], x[f"{kind}_database"], remove_alleles=True, legacy=legacy), axis=1)
    ignoring["name_similarity_without_alleles"] = ignoring.apply(
        lambda x: _levenshtein.similarity(x["cleaner_source_name"], x["cleaner_target_name"]), axis=1)
    fit_without = np.polyfit(ignoring["sequence_identity"], ignoring["name_similarity_without_alleles"], 1)
    ignoring = _fitted(ignoring, "name_similarity_without_alleles", fit_without)
    # ADR 0006: the allele-free threshold comes from its own regression; legacy used the full-name fit
    predicted_without = np.polyval(fit if legacy else fit_without, IDENTITY_THRESHOLD)
    high_without = ignoring.query("sequence_identity >= @IDENTITY_THRESHOLD")
    low_without = ignoring.query(
        "sequence_identity >= @IDENTITY_THRESHOLD and name_similarity_without_alleles <= @predicted_without")
    return {
        "distance": distance, "skipped_blaPAO": skipped, "linear_fit": fit.tolist(), "fitted": fitted,
        "predicted_similarity": float(predicted), "high_identity": high, "low_name_similarity": low,
        "distance_ignoring_alleles": ignoring, "linear_fit_without_alleles": fit_without.tolist(),
        "predicted_similarity_without_alleles": float(predicted_without),
        "high_identity_without_alleles": high_without, "low_name_similarity_without_alleles": low_without,
    }
