"""Real families drawn in the hero figure (spec #34, ticket #41), chosen by stated rules, not by hand.

All rules work on equivalence groups (proteins joined by reciprocal best hits, analysis.groups) of the release graph and
are deterministic: candidates are sorted and the first (or, for 'different', the lowest name agreement) wins.

  agree        one node per database, identical names                    -> the first by name
  punctuation  one node per database, names differ only in case or punctuation (equal after cleaning)
  allele       one node per database, names differ only in the allele number (equal once alleles are ignored)
  different    one node per database, names differ otherwise            -> the lowest agreement of the cleaned names
  family       the first group of exactly 6 proteins (many alleles of one family)
  one_way      a CARD+NCBI pair whose two members have a one-way hit to the same single ResFinder node and no RBH to ResFinder

For the first four, every RBH edge inside the group must have identity and coverage of at least 0.999 (identical
proteins) and every name at most MAX_NAME_LENGTH characters (so the labels fit the figure).
"""
import itertools

from analysis import groups, matches

MAX_NAME_LENGTH = 14
MIN_IDENTICAL = 0.999
ORDER = ["card", "ncbi", "resfinder"]
FAMILY_SIZE = 6


def _database(node_id):
    return node_id.split(":", 1)[0]


def _rank(node_id):
    return ORDER.index(_database(node_id))


def _name(graph, node_id):
    return graph.nodes[node_id]["name"]


def _rbh_pairs(graph, nodes):
    """RBH pairs inside a node set, once each, source in database order, sorted."""
    pairs = []
    for u, v in itertools.permutations(sorted(nodes), 2):
        data = graph.get_edge_data(u, v)
        if data and data["type"] == "RBH" and (_rank(u), u) < (_rank(v), v):
            pairs.append((u, v))
    return sorted(pairs, key=lambda e: (_rank(e[0]), _rank(e[1]), e))


def _identical(graph, pairs):
    return all(min(graph.get_edge_data(u, v)["identity"], graph.get_edge_data(u, v)["coverage"],
                   graph.get_edge_data(v, u)["identity"], graph.get_edge_data(v, u)["coverage"]) >= MIN_IDENTICAL
               for u, v in pairs)


def _names(graph, nodes, **kwargs):
    return [matches.clean_name(_name(graph, n), _database(n), **kwargs) for n in nodes]


def _agreement(cleaned):
    return min(matches._levenshtein.similarity(a, b) for a, b in itertools.combinations(cleaned, 2))


def select_families(graph):
    found = groups.equivalence_groups(graph)
    triples = []
    for group in found:
        nodes = sorted(group, key=lambda n: (_rank(n), n))
        if len(nodes) != 3 or [_database(n) for n in nodes] != ORDER:
            continue
        if any(len(_name(graph, n)) > MAX_NAME_LENGTH for n in nodes):
            continue
        pairs = _rbh_pairs(graph, nodes)
        if not pairs or not _identical(graph, pairs):
            continue
        triples.append((tuple(_name(graph, n) for n in nodes), nodes, pairs))
    triples.sort(key=lambda t: (t[0], t[1]))

    def pick(test, key=lambda t: 0):
        candidates = [t for t in triples if test(t)]
        return min(candidates, key=lambda t: (key(t), t[0], t[1])) if candidates else None

    raw = lambda t: len(set(t[0])) == 1
    cleaned = lambda t: len(set(_names(graph, t[1]))) == 1
    free = lambda t: len(set(_names(graph, t[1], remove_alleles=True))) == 1
    chosen = [
        ("agree", "Names agree", pick(raw)),
        ("punctuation", "Case and punctuation differ", pick(lambda t: not raw(t) and cleaned(t))),
        ("allele", "Allele number differs", pick(lambda t: not cleaned(t) and free(t))),
        ("different", "Different names, identical proteins", pick(lambda t: not free(t), key=lambda t: _agreement(_names(graph, t[1])))),
    ]
    result = [{"key": k, "title": title, "nodes": t[1], "edges": t[2]} for k, title, t in chosen if t]

    for group in found:
        nodes = sorted(group, key=lambda n: (_rank(n), n))
        if len(nodes) == FAMILY_SIZE and all(len(_name(graph, n)) <= MAX_NAME_LENGTH for n in nodes):
            result.append({"key": "family", "title": "Many alleles of one family", "nodes": nodes, "edges": _rbh_pairs(graph, nodes)})
            break

    for group in found:
        nodes = sorted(group, key=lambda n: (_rank(n), n))
        if len(nodes) == 2 and [_database(n) for n in nodes] == ["card", "ncbi"]:
            card, ncbi = nodes
            if any(len(_name(graph, n)) > MAX_NAME_LENGTH for n in nodes):
                continue
            targets = [{t for t in graph.successors(n) if _database(t) == "resfinder"
                        and graph.get_edge_data(n, t)["type"] == "OWH"} for n in nodes]
            if len(targets[0]) == 1 and targets[0] == targets[1]:
                (target,) = targets[0]
                if len(_name(graph, target)) <= MAX_NAME_LENGTH:
                    result.append({"key": "one_way", "title": "One-way hits", "nodes": nodes + [target],
                                   "edges": [(card, ncbi), (card, target), (ncbi, target)]})
                    break
    return result
