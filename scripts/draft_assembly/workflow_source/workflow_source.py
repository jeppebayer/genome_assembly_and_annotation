#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
from gwf.executors import Conda
import os, yaml, glob, sys
from workflow_templates import *

def draft_assembly_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Creates draft assembly PacBio HiFi data.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	CONDA_ENV_ASSEMBLY: str = CONFIG['condaEnvironment01']
	CONDA_ENV_QC: str | None = CONFIG['condaEnvironment02'] if CONFIG['condaEnvironment02'] else None
	SPECIES_NAME: str = CONFIG['speciesName']
	TAXONOMY: str = CONFIG['taxonomicGroup'].lower().replace(' ', '_') if CONFIG['taxonomicGroup'] else CONFIG['taxonomicGroup']
	WORK_DIR: str =  CONFIG['workingDirectoryPath'][:len(CONFIG['workingDirectoryPath']) - 1] if CONFIG['workingDirectoryPath'].endswith('/') else CONFIG['workingDirectoryPath']
	OUTPUT_DIR: str | None = (CONFIG['outputDirectoryPath'][:len(CONFIG['outputDirectoryPath']) - 1] if CONFIG['outputDirectoryPath'].endswith('/') else CONFIG['outputDirectoryPath']) if CONFIG['outputDirectoryPath'] else None
	HIFI_SEQ: list = CONFIG['pacbioHifiSequenceFiles']
	USE_HIC: int = CONFIG['useHicData']
	if USE_HIC == 1 or USE_HIC == 3:
		HIC_READS: list = CONFIG['hicSequenceFiles']
	BUSCO_DATASET: str = CONFIG['buscoDatasetName']
	POOL_SEQ_DATA: int = CONFIG['poolSeqData']
	PURGE_ROUNDS: int = CONFIG['nPurgeRounds']

	softwareList = ['hifiasm', 'busco', 'minimap2', 'purge_dups', 'seqkit', 'hmmsearch', 'metaeuk', 'miniprot']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT},
		executor=Conda(CONDA_ENV_ASSEMBLY)
	)
	
	topDir = f'{WORK_DIR}/{TAXONOMY}/{SPECIES_NAME.replace(" ", "_")}/draft_assembly' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/draft_assembly'
	topOut = f'{OUTPUT_DIR}/{TAXONOMY}/draft_assembly/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/draft_assembly/{SPECIES_NAME.replace(" ", "_")}'

	removeAdapters = gwf.target_from_template(
		name=f'{speciesAbbreviation(SPECIES_NAME)}_adapterremoval',
		template=hifiadapterfilt(
			pacbioHifiFiles=HIFI_SEQ,
			outputDirectory=topDir,
			speciesName = SPECIES_NAME,
			environment=CONDA_ENV_ASSEMBLY
		)
	)
	
	# Create draft assembly using only PacBio HiFi data
	if USE_HIC == 2 or USE_HIC == 3:
		hifiasmPrimaryAssembly = gwf.target_from_template(
			name=f'{speciesAbbreviation(SPECIES_NAME)}_hifiasm_primary',
			template=hifiasm_primary(
				hifiSequenceFile=removeAdapters.outputs['filt'],
				outputDirectory=topDir,
				speciesName=SPECIES_NAME,
				similarityThreshold=0.1 if POOL_SEQ_DATA == 1 else 0.75,
				purgeLevel=3 if POOL_SEQ_DATA == 1 else 1,
				environment=CONDA_ENV_ASSEMBLY
			)
		)

		buscoPrimary = gwf.target_from_template(
		name=f'{speciesAbbreviation(SPECIES_NAME)}_busco_primary',
		template=busco_genome(
			genomeAssemblyFile=hifiasmPrimaryAssembly.outputs['fasta'],
			buscoDataset=BUSCO_DATASET,
			environment=CONDA_ENV_QC if CONDA_ENV_QC else CONDA_ENV_ASSEMBLY
			)
		)

		for i in range(1, PURGE_ROUNDS + 1):
			if i == 1:
				inFile = hifiasmPrimaryAssembly.outputs['fasta']
			else:
				inFile = step3Primary.outputs['purged']
			step1Primary = gwf.target_from_template(
				name=f'{speciesAbbreviation(SPECIES_NAME)}_primary_pdups_s1_r{i:02}',
				template=purge_dups_1_map_hifi_to_genome(
					gemoneAssemblyFile=inFile,
					hifiSequenceFile=removeAdapters.outputs['filt'],
					outputDirectory=topDir,
					speciesName=SPECIES_NAME,
					directoryAddition='primary',
					roundNumber=i,
					environment=CONDA_ENV_ASSEMBLY
				)
			)
			step2Primary = gwf.target_from_template(
				name=f'{speciesAbbreviation(SPECIES_NAME)}_primary_pdups_s2_r{i:02}',
				template=purge_dups_2_map_to_self(
					genomeAssemblyFile=inFile,
					outputDirectory=topDir,
					speciesName=SPECIES_NAME,
					directoryAddition='primary',
					roundNumber=i,
					environment=CONDA_ENV_ASSEMBLY
				)
			)
			step3Primary = gwf.target_from_template(
				name=f'{speciesAbbreviation(SPECIES_NAME)}_primary_pdups_s3_r{i:02}',
				template=purge_dups_3_purge_duplicates(
					pbStatFile=step1Primary.outputs['stat'],
					pbBaseCovFile=step1Primary.outputs['cov'],
					selfAlignmentPaf=step2Primary.outputs['paf'],
					genomeAssemblyFile=inFile,
					outputDirectory=topDir,
					speciesName=SPECIES_NAME,
					directoryAddition='primary',
					roundNumber=i,
					environment=CONDA_ENV_ASSEMBLY
				)
			)
			purgeBuscoPrimary = gwf.target_from_template(
				name=f'{speciesAbbreviation(SPECIES_NAME)}_primary_pdups_busco_r{i:02}',
				template=busco_genome(
					genomeAssemblyFile=step3Primary.outputs['purged'],
					buscoDataset=BUSCO_DATASET,
					environment=CONDA_ENV_QC if CONDA_ENV_QC else CONDA_ENV_ASSEMBLY
				)
			)

	# Create draft assembly with PacBio HiFi data and Hi-C data
	if USE_HIC == 1 or USE_HIC == 3:

		hifiasmHicAssembly = gwf.target_from_template(
			name=f'{speciesAbbreviation(SPECIES_NAME)}_hifiasm_hic',
			template=hifiasm_hic(
				hifiSequenceFile=removeAdapters.outputs['filt'],
				hicSequenceFiles=HIC_READS,
				outputDirectory=topDir,
				speciesName=SPECIES_NAME,
				similarityThreshold=0.1 if POOL_SEQ_DATA == 1 else 0.75,
				purgeLevel=3 if POOL_SEQ_DATA == 1 else 1,
				environment=CONDA_ENV_ASSEMBLY
			)
		)

		buscoHic = gwf.target_from_template(
			name=f'{speciesAbbreviation(SPECIES_NAME)}_busco_hic',
			template=busco_genome(
				genomeAssemblyFile=hifiasmHicAssembly.outputs['fasta'][0],
				buscoDataset=BUSCO_DATASET,
				environment=CONDA_ENV_QC if CONDA_ENV_QC else CONDA_ENV_ASSEMBLY
			)
		)

		for i in range(1, PURGE_ROUNDS + 1):
			if i == 1:
				inFile = hifiasmHicAssembly.outputs['fasta'][0]
			else:
				inFile = step3Hic.outputs['purged']
			step1Hic = gwf.target_from_template(
				name=f'{speciesAbbreviation(SPECIES_NAME)}_hic_pdups_s1_r{i:02}',
				template=purge_dups_1_map_hifi_to_genome(
					gemoneAssemblyFile=inFile,
					hifiSequenceFile=removeAdapters.outputs['filt'],
					outputDirectory=topDir,
					speciesName=SPECIES_NAME,
					directoryAddition='hic',
					roundNumber=i,
					environment=CONDA_ENV_ASSEMBLY
				)
			)
			step2Hic = gwf.target_from_template(
				name=f'{speciesAbbreviation(SPECIES_NAME)}_hic_pdups_s2_r{i:02}',
				template=purge_dups_2_map_to_self(
					genomeAssemblyFile=inFile,
					outputDirectory=topDir,
					speciesName=SPECIES_NAME,
					directoryAddition='hic',
					roundNumber=i,
					environment=CONDA_ENV_ASSEMBLY
				)
			)
			step3Hic = gwf.target_from_template(
				name=f'{speciesAbbreviation(SPECIES_NAME)}_hic_pdups_s3_r{i:02}',
				template=purge_dups_3_purge_duplicates(
					pbStatFile=step1Hic.outputs['stat'],
					pbBaseCovFile=step1Hic.outputs['cov'],
					selfAlignmentPaf=step2Hic.outputs['paf'],
					genomeAssemblyFile=inFile,
					outputDirectory=topDir,
					speciesName=SPECIES_NAME,
					directoryAddition='hic',
					roundNumber=i,
					environment=CONDA_ENV_ASSEMBLY
				)
			)
			purgeBusco = gwf.target_from_template(
				name=f'{speciesAbbreviation(SPECIES_NAME)}_hic_pdups_busco_r{i:02}',
				template=busco_genome(
					genomeAssemblyFile=step3Hic.outputs['purged'],
					buscoDataset=BUSCO_DATASET,
					environment=CONDA_ENV_QC if CONDA_ENV_QC else CONDA_ENV_ASSEMBLY
				)
			)

	# if ASSEMBLY:
	#     merqury_qc = gwf.target_from_template(
	#         name=f'{speciesAbbreviation(SPECIES_NAME)}_merqury',
	#         template=merqury(
	#             genome_assembly_file=ASSEMBLY,
	#             pacbio_hifi_reads=remove_adapters.outputs['filt'],
	#             output_directory=top_dir
	#         )
	#     )

	#     blobtools_blastn_search = gwf.target_from_template(
	#         name=f'{speciesAbbreviation(SPECIES_NAME)}_blobtools_blastn',
	#         template=blobtools_blastn(
	#             genome_assembly_file=ASSEMBLY
	#         )
	#     )

	#     blobtools_diamond_search = gwf.target_from_template(
	#         name=f'{speciesAbbreviation(SPECIES_NAME)}_blobtools_diamond',
	#         template=blobtools_diamond(
	#             genome_assembly_file=ASSEMBLY
	#         )
	#     )

	#     blobtools_coverage_alignment = gwf.target_from_template(
	#         name=f'{speciesAbbreviation(SPECIES_NAME)}_blobtools_coverage',
	#         template=blobtools_coverage(
	#             genome_assembly_file=ASSEMBLY,
	#             pacbio_hifi_reads=remove_adapters.outputs['filt']
	#         )
	#     )

	#     blobtools = gwf.target_from_template(
	#         name=f'{speciesAbbreviation(SPECIES_NAME)}_blobtools_blobdir',
	#         template=blobtools_blobdir(
	#             genome_assembly_file=ASSEMBLY,
	#             species_name=SPECIES_NAME,
	#             blastn_result_file=blobtools_blastn_search.outputs['blast'],
	#             diamond_result_file=blobtools_diamond_search.outputs['diamond'],
	#             coverage_alignment_file=blobtools_coverage_alignment.outputs['alignment'],
	#             busco_full_table_file="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Isotoma_viridis/draft_assembly/purge_dups_primary/04/busco_IsoVir.purged.fa/run_arthropoda_odb10/full_table.tsv"
	#         )
	#     )

	with open(f'softwareVersions.tsv', 'w') as outfile:
		outfile.write(software_versions_to_string(software_versions([CONDA_ENV_ASSEMBLY, CONDA_ENV_QC], softwareList)))

	print(f'Intermediary files will be place at: {topDir}/')
	print(f'Output files will be placed at: {topOut if OUTPUT_DIR else topDir}/')

	return gwf