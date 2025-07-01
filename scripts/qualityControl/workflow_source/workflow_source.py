#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
from gwf.executors import Conda
import os, yaml, glob, sys
from workflow_templates import *
import importlib.util

spec = importlib.util.spec_from_file_location('info', f'{os.path.dirname(os.path.realpath(__file__))}/software/taxon_info/taxonInfo.py')
info = importlib.util.module_from_spec(spec)
spec.loader.exec_module(info)

def assembly_quality_control_workflow(configFile: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: description
	
	:param str configFile:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(configFile))
	ACCOUNT: str = CONFIG['account']
	CONDA_ENV: str = CONFIG['condaEnvironment01']
	TAXONOMY: str | None = CONFIG['taxonomicGroup'].lower().replace(' ', '_') if CONFIG['taxonomicGroup'] else None
	SPECIES_NAME: str = CONFIG['speciesName']
	WORK_DIR: str =  CONFIG['workingDirectoryPath'][:len(CONFIG['workingDirectoryPath']) - 1] if CONFIG['workingDirectoryPath'].endswith('/') else CONFIG['workingDirectoryPath']
	OUTPUT_DIR: str | None = (CONFIG['outputDirectoryPath'][:len(CONFIG['outputDirectoryPath']) - 1] if CONFIG['outputDirectoryPath'].endswith('/') else CONFIG['outputDirectoryPath']) if CONFIG['outputDirectoryPath'] else None
	ASSEMBLY_FILE: str = CONFIG['genomeAssemblyFile']
	SAMPLE_SETUP: dict = CONFIG['sampleSetup']
	SAMPLE_NAME: str = SAMPLE_SETUP['sampleName']
	DATATYPE: str = SAMPLE_SETUP['datatype']
	SEQUENCING_FILES: list = SAMPLE_SETUP['sequencingFiles']
	DATABASES: dict = CONFIG['databases']
	BUSCO_PATH: str = DATABASES['busco']['path']
	REQUESTED_BUSCO_LINEAGES: str | None = ','.join(DATABASES['busco']['requestedLineages']) if DATABASES['busco']['requestedLineages'][0] else None
	UNIPROT: str = DATABASES['uniprotPath']
	NCBI_TAXDUMP: str = DATABASES['ncbiTaxdumpPath']
	NCBI_NT: str = DATABASES['ncbiNtPath']
	
	BASAL_LINEAGES = ['eukaryota_odb12', 'bacteria_odb12', 'archaea_odb12']

	softwareList = ['fasta_windows', 'busco', 'blobtoolkit', 'blast', 'diamond', 'minimap2', 'samtools']

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------

	if not os.path.exists(f'./{speciesAbbreviation(SPECIES_NAME)}.info.yml'):
		print(f'Creating {speciesAbbreviation(SPECIES_NAME)}.info.yml...')
		taxonInfo = info.get_taxon_info(SPECIES_NAME, NCBI_TAXDUMP)
		classification = info.get_classification(taxonInfo)
		odbList = info.get_odb_list(taxonInfo, REQUESTED_BUSCO_LINEAGES)
		assemblyInfo = info.get_assembly_info(ASSEMBLY_FILE)
		info.print_yaml(speciesAbbreviation(SPECIES_NAME), taxonInfo, assemblyInfo, classification, odbList)
		print('File has been created... Running workflow...')

	INFO = yaml.safe_load(open(f'./{speciesAbbreviation(SPECIES_NAME)}.info.yml'))
	BUSCO_LINEAGES = INFO['buscoLineages']
	TAX_ID = INFO['taxonomy']['taxid']

	orderedLineages = {index: {'odb': lineage, 'fullTable': str} for index, lineage in enumerate(BUSCO_LINEAGES)}
	lastIndex = max(list(orderedLineages.keys()))
	for basal_lineage in BASAL_LINEAGES:
		if basal_lineage not in BUSCO_LINEAGES:
			lastIndex += 1
			orderedLineages[lastIndex] = {'odb': basal_lineage, 'fullTable': str}

	gwf = Workflow(
		defaults={'account': ACCOUNT},
        executor=Conda(CONDA_ENV)
	)
	
	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/assemblyQC/{os.path.basename(ASSEMBLY_FILE)}' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/assemblyQC/{os.path.basename(ASSEMBLY_FILE)}'
	topOut = f'{OUTPUT_DIR}/assemblyQC/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/{os.path.basename(ASSEMBLY_FILE)}' if TAXONOMY else f'{OUTPUT_DIR}/assemblyQC/{SPECIES_NAME.replace(" ", "_")}/{os.path.basename(ASSEMBLY_FILE)}'

	# genomeSize = gwf.target_from_template(
		# name=f'genome_size',
		# template=genome_size(
			# genomeAssemblyFile=ASSEMBLY_FILE,
			# outputDirectory=topDir,
			# environment=CONDA_ENV
		# )
	# )

	fastaWindows = gwf.target_from_template(
		name=f'fasta_windows',
		template=fasta_windows(
			genomeAssemblyFile=ASSEMBLY_FILE,
			outputDirectory=topDir,
			environment=CONDA_ENV
		)
	)

	blobtoolkitAlignment = gwf.target_from_template(
		name=f'alignment_hifi_to_assembly',
		template=blobtoolkit_alignment(
			genomeAssemblyFile=ASSEMBLY_FILE,
			pacbioHifiReads=SEQUENCING_FILES,
			outputDirectory=topDir,
			environment=CONDA_ENV
		)
	)

	samtoolsFlagstat = gwf.target_from_template(
		name=f'samtools_flagstat',
		template=samtools_flagstat(
			alignmnetBamFile=blobtoolkitAlignment.outputs['alignment'],
			outputDirectory=topDir,
			environment=CONDA_ENV
		)
	)

	blobtoolkitCoverage = gwf.target_from_template(
		name=f'blobtoolkit_coverage',
		template=blobtoolkit_coverage(
			hifiToAssemblyBam=blobtoolkitAlignment.outputs['alignment'],
			outputDirectory=topDir,
			environment=CONDA_ENV
		)
	)

	for index, lineage in enumerate([rank['odb'] for rank in orderedLineages.values()]):
		buscoGenome = gwf.target_from_template(
			name=f'busco_genome_{lineage}',
			template=busco_genome(
				genomeAssemblyFile=ASSEMBLY_FILE,
				buscoLineage=lineage,
				buscoDownloadPath=BUSCO_PATH,
				outputDirectory=topDir,
				environment=CONDA_ENV
			)
		)

		orderedLineages[index]['fullTable'] = buscoGenome.outputs['fulltable']

		if lineage in BASAL_LINEAGES:
			blobtoolkitExtractBuscoGenes = gwf.target_from_template(
				name=f'extract_busco_genes_{lineage}',
				template=blobtoolkit_extract_busco_genes(
					buscoFullTableTsv=buscoGenome.outputs['fulltable'],
					outputPrefix=f'{os.path.basename(os.path.splitext(os.path.splitext(ASSEMBLY_FILE)[0])[0]) if ASSEMBLY_FILE.endswith(".gz") else os.path.basename(os.path.splitext(ASSEMBLY_FILE)[0])}.{lineage}',
					outputDirectory=topDir,
					environment=CONDA_ENV
				)
			)

			diamondBlastp = gwf.target_from_template(
				name=f'diamond_blastp_{lineage}',
				template=diamond_blastp(
					queryFileFasta=blobtoolkitExtractBuscoGenes.outputs['fasta'],
					excludeTaxon=TAX_ID,
					outputDirectory=topDir,
					environment=CONDA_ENV,
					diamondDatabaseFile=UNIPROT
				)
			)

	blobtoolkitCountBuscoGenes = gwf.target_from_template(
		name=f'count_busco_genes',
		template=blobtoolkit_count_busco_genes(
			buscoTablesFull=[rank['fullTable'] for rank in orderedLineages.values()],
			filename=os.path.basename(os.path.splitext(os.path.splitext(ASSEMBLY_FILE)[0])[0]) if ASSEMBLY_FILE.endswith('.gz') else os.path.basename(os.path.splitext(ASSEMBLY_FILE)[0]),
			windowsBedFile=fastaWindows.outputs['bed'],
			outputDirectory=topDir,
			environment=CONDA_ENV
		)
	)

	blobtoolkitCollateStats = gwf.target_from_template(
		name=f'collate_stats',
		template=blobtoolkit_collate_stats(
			countBuscoGenesFile=blobtoolkitCountBuscoGenes.outputs['count'],
			windowFreqFile=fastaWindows.outputs['freq'],
			windowMononucFile=fastaWindows.outputs['mono'],
			regionsCoverageFile=blobtoolkitCoverage.outputs['bed'],
			filename=os.path.basename(os.path.splitext(os.path.splitext(ASSEMBLY_FILE)[0])[0]) if ASSEMBLY_FILE.endswith('.gz') else os.path.basename(os.path.splitext(ASSEMBLY_FILE)[0]),
			outputDirectory=topDir,
			environment=CONDA_ENV
		)
	)

	diamondBlastx = gwf.target_from_template(
		name=f'diamond_blastx',
		template=diamond_blastx(
			queryFileFasta=ASSEMBLY_FILE,
			buscoTableFull=orderedLineages[0]['fullTable'],
			buscoLineage=orderedLineages[0]["odb"],
			excludeTaxon=TAX_ID,
			outputDirectory=topDir,
			environment=CONDA_ENV,
			diamondDatabaseFile=UNIPROT
		)
	)

	diamondBlastxNoHits = gwf.target_from_template(
		name=f'diamond_blastx_no_hits',
		template=diamond_blastx_no_hits(
			queryFileFasta=ASSEMBLY_FILE,
			blastxResults=diamondBlastx.outputs['blast'],
			outputDirectory=topDir,
			environment=CONDA_ENV
		)
	)

	NcbiBlastnWithoutTaxon = gwf.target_from_template(
		name=f'ncbi_blastn_without_taxon_{TAX_ID}',
		template=ncbi_blastn(
			queryFileFasta=ASSEMBLY_FILE,
			outputDirectory=topDir,
			excludeTaxon=TAX_ID,
			environment=CONDA_ENV,
			ncbiBlastDatabase=NCBI_NT
		)
	)

	ncbiBlastnWithTaxon = gwf.target_from_template(
		name=f'ncbi_blastn_with_taxon_{TAX_ID}',
		template=ncbi_blastn(
			queryFileFasta=ASSEMBLY_FILE,
			outputDirectory=topDir,
			excludeTaxon='',
			environment=CONDA_ENV,
			ncbiBlastDatabase=NCBI_NT
		)
	)
	
	with open(f'softwareVersions.tsv', 'w') as outfile:
		outfile.write(software_versions_to_string(software_versions([CONDA_ENV], softwareList)))

	print(f'Intermediary files will be place at: {topDir}/')
	print(f'Output files will be placed at: {topOut if OUTPUT_DIR else topDir}/')

	return gwf