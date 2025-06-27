#!/bin/env python3

querySpecies = 'Lepidocyrtus lignorum'
taxdumpPath = '/faststorage/project/EcoGenetics/databases/NCBI_Taxdump_062025'

nodes = taxdumpPath + '/' + 'nodes.dmp'
fullnamelineage = taxdumpPath + '/' + 'fullnamelineage.dmp'
rankedlineage = taxdumpPath + '/' + 'rankedlineage.dmp'
taxidlineage = taxdumpPath + '/' + 'taxidlineage.dmp'

odbMapping = {'odb10': '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/scripts/qualityControl/workflow_source/software/mapping_taxids-busco_dataset_name.odb10.2019-12-16.txt',
			  'odb12': '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/scripts/qualityControl/workflow_source/software/mapping_taxids-busco_dataset_name.odb12.2025-01-15.txt'}

RANKS = [
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "domain",
]

taxonInfo = {'taxonId': int, 'name': str, 'rank': str, 'lineage': list}

with open(fullnamelineage, 'r') as infile:
	for line in infile:
		entry = line.replace('\t', '').strip().split('|')[:3]
		if entry[1] == querySpecies:
			taxonInfo['taxonId'] = entry[0]
			taxonInfo['name'] = entry[1]
			taxonInfo['lineage'] = [{'taxonId': int, 'name': lineage.strip(), 'rank': str} for lineage in reversed(entry[2].split(';')[1:-1])]
			break

with open(taxidlineage, 'r') as infile:
	for line in infile:
		entry = line.replace('\t', '').strip().split('|')[:2]
		if entry[0] == taxonInfo['taxonId']:
			ancestors = entry[1].strip().split(' ')[1:]
			for index, id in enumerate(reversed(ancestors)):
				taxonInfo['lineage'][index]['taxonId'] = id
			break

with open(nodes, 'r') as infile:
	for line in infile:
		entry = line.replace('\t', '').strip().split('|')[:3]
		if entry[0] == taxonInfo['taxonId']:
			taxonInfo['rank'] = entry[2].upper()
			break

with open(nodes, 'r') as infile:
	for line in infile:
		entry = line.replace('\t', '').strip().split('|')[:3]
		for index, id in enumerate(taxonInfo['lineage']):
			if entry[0] == id['taxonId']:
				taxonInfo['lineage'][index]['rank'] = entry[2].upper()
				break

classification = {ancestor['rank'].lower(): ancestor['name'] for ancestor in taxonInfo['lineage'] if ancestor['rank'].lower() in RANKS}


print(taxonInfo)
print(classification)

# taxon:
#   class: Collembola
#   family: Entomobryidae
#   genus: Lepidocyrtus
#   kingdom: Metazoa
#   name: Lepidocyrtus lignorum
#   order: Entomobryomorpha
#   phylum: Arthropoda
#   superkingdom: Eukaryota
#   taxid: '707889'