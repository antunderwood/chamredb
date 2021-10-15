import pronto
import json
import os
import re
import argparse

from pronto import relationship

parser = argparse.ArgumentParser(description='Process input/output files')
parser.add_argument('--input', type=str, help='input file name')
parser.add_argument('--output', type=str, help='output file name')

args = parser.parse_args()

if args.input:
	metadata_file = args.input
else:
	metadata_file =  os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"..", "chamredb", "data", "db_metadata",
		"card.metadata.obo"
	)

if args.output:
	out_path = args.output
else:
	out_path =  os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"..", "chamredb", "data", "db_metadata",
		"card.metadata.json"
	)

o = pronto.Ontology(metadata_file)

metadata = {}
for term in o.terms():
    metadata[term.id] = {'name': term.name}
    # add is_a superclass
    if len(list(term.superclasses(distance=1, with_self=False))) > 0:
        is_a = ','.join([f"{term.id}:{term.name}" for term in term.superclasses(distance=1, with_self=False)])
        metadata[term.id]['is_a'] = is_a
    # parse relationship
    relationships = term.relationships
    for relationship in relationships:
        if relationship.name.startswith('confers_resistance'):
            antibiotics = set([term.name for term in relationships[relationship]])
            if "phenotype" not in metadata[term.id]:
                metadata[term.id]["phenotype"] = {}
            metadata[term.id]["phenotype"][f"{relationship.name.replace('_', ' ')}"] = antibiotics

    # add publication metadata
    if len(list(term.definition.xrefs)) > 0:
        metadata[term.id]["PMID"] = []
        for item in term.definition.xrefs:
            metadata[term.id]["PMID"].append(item.id.replace("PMID:", ""))

# second pass to pull out extra confers resistance from is_a relationships
for id in metadata:
    if "is_a" in metadata[id]:
        is_a_ARO_id = ":".join(metadata[id]["is_a"].split(":")[0:2])
        if "phenotype" in metadata[is_a_ARO_id]:
            if "phenotype" not in metadata[id]:
                metadata[id]["phenotype"] = {}
            for phenotype_type in metadata[is_a_ARO_id]["phenotype"]:
                if phenotype_type in metadata[id]["phenotype"]:
                    metadata[id]["phenotype"][phenotype_type].update(metadata[is_a_ARO_id]["phenotype"][phenotype_type])
                else:
                    metadata[id]["phenotype"][phenotype_type] = metadata[is_a_ARO_id]["phenotype"][phenotype_type]

# finally convert phenotypes to strings
for id in metadata:
    if "phenotype" in  metadata[id]:
        for phenotype_type in metadata[id]["phenotype"]:
            if phenotype_type not in ["confers resistance to antibiotic", "confers resistance to drug class"]:
                print(phenotype_type)
            metadata[id]["phenotype"][phenotype_type] = ",".join(sorted(list(metadata[id]["phenotype"][phenotype_type])))
        primary_type = "confers resistance to antibiotic"
        secondary_type = "confers resistance to drug class"
        if primary_type in metadata[id]["phenotype"]:
            if secondary_type in metadata[id]["phenotype"]:
                metadata[id]["additional_phenotype"] = f"{secondary_type}: {metadata[id]['phenotype'][secondary_type]}"
            metadata[id]["phenotype"] = f"{primary_type}: {metadata[id]['phenotype'][primary_type]}"
        else:
            metadata[id]["phenotype"] = f"{secondary_type}: {metadata[id]['phenotype'][secondary_type]}"
            

with open(out_path, "w") as out_file:
    out_file.write(
        json.dumps(metadata, sort_keys=True, indent=2)
    )
