#!/bin/env python3
import sys, os, yaml, gzip, argparse

RANKS = [
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "superkingdom",
]

def parse_args(args = None):
	description = 'Create configuration file blobtoolkit analysis.'

	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--query', help='Scientific name of query species', required=True)
	parser.add_argument('--assembly', help='FASTA format genome assembly file', required=True)
	parser.add_argument('--taxdump', help='Path to NCBI Taxdump directory', required=True)
	parser.add_argument('--uniprot', help='Path to UniProt reference proteomes database', required=True)
	parser.add_argument('--nt', help='Path to NCBI nt database', required=True)
	parser.add_argument('--busco', help='Path to BUSCO database', required=True)
	parser.add_argument('--lineage', help='Comma-separated list of requested BUSCO lineages', default=None)
	parser.add_argument('--readlayout', action='append', choices=['single', 'paired'], type=str, help='Sequencing layout', required=True)
	parser.add_argument('--readtype', action='append', help='Type of read set', required=True)
	parser.add_argument('--reads', action='append', help='Comma-separated path to reads file set', required=True)
	parser.add_argument('--prefix', help='Output prefix', required=True)
	args = parser.parse_args(args)
	if not os.path.exists(args.taxdump):
		print('TAXDUMP must be an existing directory', file=sys.stderr)
		sys.exit(1)
	if not os.path.exists(args.taxdump):
		print('UNIPROT must be an existing directory', file=sys.stderr)
		sys.exit(1)
	if not os.path.exists(args.taxdump):
		print('NT must be an existing directory', file=sys.stderr)
		sys.exit(1)
	if not os.path.exists(args.taxdump):
		print('BUSCO must be an existing directory', file=sys.stderr)
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
	taxonInfo = {'taxid': int, 'name': str, 'rank': str, 'lineage': list}
	with open(fullnamelineage, 'r') as infile:
		for line in infile:
			entry = line.replace('\t', '').strip().split('|')[:3]
			if entry[1] == querySpecies:
				taxonInfo['taxid'] = entry[0]
				taxonInfo['name'] = entry[1]
				taxonInfo['lineage'] = [{'taxid': int, 'name': lineage.strip(), 'rank': str} for lineage in reversed(entry[2].split(';')[1:-1])]
				break
	if taxonInfo['taxid'] == int:
		print(f'Cannot find taxon information for \'{querySpecies}\'... Trying \'unclassified {querySpecies.split()[0]}\'...', file=sys.stderr)
		with open(fullnamelineage, 'r') as infile:
			for line in infile:
				entry = line.replace('\t', '').strip().split('|')[:3]
				if entry[1] == f'unclassified {querySpecies.split()[0]}':
					taxonInfo['taxid'] = entry[0]
					taxonInfo['name'] = f'{querySpecies} ({entry[1]})'
					taxonInfo['lineage'] = [{'taxid': int, 'name': lineage.strip(), 'rank': str} for lineage in reversed(entry[2].split(';')[1:-1])]
					break
		if taxonInfo['taxid'] == int:
			print(f'ERROR: Was not able to retrieve taxonomic information for \'{querySpecies}\'', file=sys.stderr)
			sys.exit(2)

	with open(taxidlineage, 'r') as infile:
		for line in infile:
			entry = line.replace('\t', '').strip().split('|')[:2]
			if entry[0] == taxonInfo['taxid']:
				ancestors = entry[1].strip().split(' ')[1:]
				for index, id in enumerate(reversed(ancestors)):
					taxonInfo['lineage'][index]['taxid'] = id
				break

	with open(nodes, 'r') as infile:
		for line in infile:
			entry = line.replace('\t', '').strip().split('|')[:3]
			if entry[0] == taxonInfo['taxid']:
				taxonInfo['rank'] = entry[2].upper()
				break

	with open(nodes, 'r') as infile:
		for line in infile:
			entry = line.replace('\t', '').strip().split('|')[:3]
			for index, id in enumerate(taxonInfo['lineage']):
				if entry[0] == id['taxid']:
					taxonInfo['lineage'][index]['rank'] = entry[2].upper()
					break
	return taxonInfo

def get_classification(taxonInfo: dict) -> dict:
	classification = {ancestor['rank'].lower(): ancestor['name'] for ancestor in taxonInfo['lineage']}
	if 'superkingdom' not in classification:
		classification['superkingdom'] = classification['domain']
	return {rank: classification[rank] for rank in RANKS if rank in RANKS}

def get_odb_list(taxonInfo: dict,
				 requestedBuscos: str,
				 odbMapping: dict = {'odb10': f'{os.path.dirname(os.path.realpath(__file__))}/taxidMapping/mapping_taxids-busco_dataset_name.odb10.2019-12-16.txt',
					   				 'odb12': f'{os.path.dirname(os.path.realpath(__file__))}/taxidMapping/mapping_taxids-busco_dataset_name.odb12.2025-01-15.txt'}
				) -> list:
	odbDict = {}
	for odb in odbMapping:
		with open(odbMapping[odb], 'r') as infile:
			odbDict[odb] = {}
			for line in infile:
				entry = line.strip().split('\t')
				odbDict[odb][entry[0]] = entry[1] + '_' + odb
	odbList = [odbDict['odb12'][ancestralTaxonId['taxid']] for ancestralTaxonId in taxonInfo['lineage'] if ancestralTaxonId['taxid'] in odbDict['odb12']]
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
	assemblyInfo = {'scaffold-count': 0, 'span': 0, 'file_size(MB)': round(os.path.getsize(assemblyFile) / 1000000), 'file': assemblyFile}
	with gzip.open(assemblyFile, 'rt') if assemblyFile.endswith('.gz') else open(assemblyFile, 'r') as infile:
		for line in infile:
			if line.strip().startswith('>'):
				assemblyInfo['scaffold-count'] += 1
			else:
				assemblyInfo['span'] += len(line.strip())
	return assemblyInfo

def base_count(readFiles: list):
	count = 0
	for file in readFiles:
		with gzip.open(file, 'rt') if file.endswith('.gz') else open(file, 'r') as infile:
			for line in infile:
				if not line.strip().startswith('>'):
					count += len(line.strip())
	return count

def get_platform(datatype: str) -> str:
    if datatype in ['ont', 'nanopore']:
        return 'OXFORD_NANOPORE'
    elif datatype.startswith('pacbio') or datatype == 'hifi':
        return 'PACBIO_SMRT'
    elif datatype in ['hic', 'illumina']:
        return 'ILLUMINA'
    else:
        return 'OTHER'

def print_yaml(prefix: str, taxonInfo: dict, assemblyInfo: dict, classification: dict, odbList: list, busco: str, readset: zip, nt: str, uniprot: str, taxdump: str):
	assemblyInfo['accession'] = 'draft'
	assemblyInfo['level'] = 'scaffold'
	data = {
		'assembly': assemblyInfo,
		'reads': {
			'paired': [],
			'single': []},
		'revision': 1,
		'taxon': {
			'taxid': taxonInfo['taxid'],
			'name': taxonInfo['name'],
			**classification
		},
		'busco': {
			'lineages': odbList,
			'lineage_dir': busco + '/lineages'},
		'similarity': {
			'blastn': {
				'name': 'nt',
				'path': nt,
				'source': 'ncbi'
			},
			'diamond_blastp': {
				'name': 'reference_proteomes',
				'path': uniprot,
				'source': 'uniprot',
				'import_max_target_seqs': 100000,
				'taxrule': 'blastp=buscogenes'
			},
			'diamond_blastx': {
				'name': 'reference_proteomes',
				'path': uniprot,
				'source': 'uniprot'
			},
			'defaults': {
				'evalue': 1e-10,
				'import_evalue': 1e-25,
				'max_target_seq': 10,
				'taxrule': 'buscogenes'
			}
		},
		'settings': {
			'blast_chunk': 100000,
			'blast_max_chunks': 10,
			'blast_min_length': 1000,
			'blast_overlap': 0,
			'stats_chunk': 1000,
			'stats_windows': [0.1, 0.01, 1, 100000, 1000000],
			'taxdump': taxdump
		},
		'version': 1
	}

	for readlayout, readtype, reads in readset:
		# readsList = reads.split(',')
		# data['reads']['paired' if len(readsList) > 1 else 'single'].append(
		# 		{
		# 			'file': reads.replace(',', ';'),
		# 			'platform': get_platform(readtype.lower())
		# 		}
		# 	)
		data['reads'][readlayout].append(
				{
					'file': reads.replace(',', ';'),
					'platform': get_platform(readtype.lower())
				}
			)

	if len(os.path.dirname(prefix)) > 0:
		os.makedirs(os.path.dirname(prefix), exist_ok=True)
	with open(f'{prefix}.info.yaml', 'w') as outfile:
		yaml.dump(data, outfile)

def main(args = None):
	args = parse_args(args)
	taxonInfo = get_taxon_info(args.query, args.taxdump)
	classification = get_classification(taxonInfo)
	odbList = get_odb_list(taxonInfo, args.lineage)
	assemblyInfo = get_assembly_info(args.assembly)
	readset = zip(args.readlayout, args.readtype, args.reads) if args.reads else []
	print_yaml(args.prefix, taxonInfo, assemblyInfo, classification, odbList, args.busco, readset, args.nt, args.uniprot, args.taxdump)

if __name__ == '__main__':
	sys.exit(main())