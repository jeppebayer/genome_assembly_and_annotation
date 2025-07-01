#!/bin/env python3
import sys, os, yaml, gzip, argparse

RANKS = [
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "domain",
]

def parse_args(args = None):
	description = 'Combine files to create input file for blobtoolkit window-stats.'

	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--query', help='Scientific name of query species', required=True)
	parser.add_argument('--taxdump', help='Path to NCBI Taxdump directory', required=True)
	parser.add_argument('--busco', help='Comma-separated list of requested BUSCO lineages', default=None)
	parser.add_argument('--assembly', help='Genome assembly file', required=True)
	parser.add_argument('--prefix', help='Output prefix', required=True)
	args = parser.parse_args(args)
	if not os.path.exists(args.taxdump):
		print('TAXDUMP must be an existing directory', file=sys.stderr)
		sys.exit(1)
	if not os.path.exists(args.assembly):
		print('ASSEMBLY must be an existing file', file=sys.stderr)
		sys.exit(1)
	return args

def get_taxon_info(querySpecies: str, taxdumpPath: str) -> dict: 
	nodes = taxdumpPath + '/' + 'nodes.dmp'
	fullnamelineage = taxdumpPath + '/' + 'fullnamelineage.dmp'
	rankedlineage = taxdumpPath + '/' + 'rankedlineage.dmp'
	taxidlineage = taxdumpPath + '/' + 'taxidlineage.dmp'
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
	return taxonInfo

def get_classification(taxonInfo: dict) -> dict:
	classification = {ancestor['rank'].lower(): ancestor['name'] for ancestor in taxonInfo['lineage'] if ancestor['rank'].lower() in RANKS}
	return classification

def get_odb_list(taxonInfo: dict,
				 requestedBuscos: str,
				 odbMapping: dict = {'odb10': f'{os.path.dirname(os.path.realpath(__file__))}/mapping_taxids-busco_dataset_name.odb10.2019-12-16.txt',
					   				 'odb12': f'{os.path.dirname(os.path.realpath(__file__))}/mapping_taxids-busco_dataset_name.odb12.2025-01-15.txt'}
				) -> list:
	odbDict = {}
	for odb in odbMapping:
		with open(odbMapping[odb], 'r') as infile:
			odbDict[odb] = {}
			for line in infile:
				entry = line.strip().split('\t')
				odbDict[odb][entry[0]] = entry[1] + '_' + odb
	odbList = [odbDict['odb12'][ancestralTaxonId['taxonId']] for ancestralTaxonId in taxonInfo['lineage'] if ancestralTaxonId['taxonId'] in odbDict['odb12']]
	validOdbs = []
	for odbSet in odbDict:
		validOdbs = validOdbs + [odb for odb in odbDict[odbSet].values()]
	validOdbs = set(validOdbs)
	if requestedBuscos:
		requestedBuscosList = requestedBuscos.split(',')
		for odb in requestedBuscosList:
			if odb not in validOdbs:
				print(f'Invalid requested BUSCO lineage: {odb}', file=sys.stderr)
				sys.exit(1)
			if odb not in odbList:
				odbList.append(odb)
	return odbList

def get_assembly_info(assemblyFile: str):
	assemblyInfo = {'nScaffolds': 0, 'genomeSize(Mb)': 0, 'fileSize(MB)': round(os.path.getsize(assemblyFile) / 1000000), 'assemblyFile': assemblyFile}
	with gzip.open(assemblyFile, 'rb') if assemblyFile.endswith('.gz') else open(assemblyFile, 'r') as infile:
		for line in infile:
			if line.strip().startswith('>'):
				assemblyInfo['nScaffolds'] += 1
			else:
				assemblyInfo['genomeSize(Mb)'] += len(line.strip())
	assemblyInfo['genomeSize(Mb)'] = round(assemblyInfo['genomeSize(Mb)'] / 1000000)
	return assemblyInfo


def print_yaml(prefix: str, taxonInfo: dict, assemblyInfo: dict, classification: dict, odbList: list):
	data = {
		'assembly': assemblyInfo,
		'taxonomy': {
			'taxid': taxonInfo['taxonId'],
			'name': taxonInfo['name'],
			**classification
		},
		'buscoLineages': odbList
	}
	if len(os.path.dirname(prefix)) > 0:
		os.makedirs(os.path.dirname(prefix), exist_ok=True)
	with open(f'{prefix}.info.yml', 'w') as outfile:
		yaml.dump(data, outfile)

def main(args = None):
	args = parse_args(args)
	taxonInfo = get_taxon_info(args.query, args.taxdump)
	classification = get_classification(taxonInfo)
	odbList = get_odb_list(taxonInfo, args.busco)
	assemblyInfo = get_assembly_info(args.assembly)
	print_yaml(args.prefix, taxonInfo, assemblyInfo, classification, odbList)

if __name__ == '__main__':
	sys.exit(main())