"""Equivalence groups: proteins joined by reciprocal best hits (spec #34, ticket #49).

A group is a connected component of the graph restricted to RBH edges; one-way hits never join groups and a node with no
RBH is a group of its own. The database profile of a group is the sorted set of databases of its nodes, for example
"card+ncbi".
"""
from collections import Counter

import networkx as nx


def _database(node_id):
    return node_id.split(":", 1)[0]


def equivalence_groups(graph):
    """List of frozensets of node ids, ordered by the smallest node id of each group (deterministic)."""
    rbh = nx.Graph()
    rbh.add_nodes_from(graph.nodes)
    rbh.add_edges_from((u, v) for u, v, d in graph.edges(data=True) if d["type"] == "RBH")
    return sorted((frozenset(c) for c in nx.connected_components(rbh)), key=min)


def profiles(found):
    """Counts of groups and of proteins by database profile, the group size distribution and the largest group."""
    group_counts, protein_counts, sizes = Counter(), Counter(), Counter()
    for group in found:
        profile = "+".join(sorted({_database(n) for n in group}))
        group_counts[profile] += 1
        protein_counts[profile] += len(group)
        sizes[len(group)] += 1
    largest = max(found, key=lambda g: (len(g), sorted(g)))
    return {
        "groups": dict(sorted(group_counts.items())), "proteins": dict(sorted(protein_counts.items())),
        "size_distribution": dict(sorted(sizes.items())), "largest_group_size": len(largest),
        "largest_group_by_database": dict(sorted(Counter(_database(n) for n in largest).items())),
    }
