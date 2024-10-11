#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def repeat_masking_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Uses RepeatModeler, RepeatMasker and data from available databases to identify and soft-mask repetitive region in a genome assembly.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	SPECIES_NAME: str = CONFIG['species_name']
	GENOME_ASSEMBLY: str = CONFIG['genome_assembly']
	REPEAT_DATABASE: str = CONFIG['repeat_database']
	WORKING_DIR: str = CONFIG['working_directory_path']
	OUTPUT_DIR: str = CONFIG['output_directory_path']
	ALT_NAME: str = CONFIG['alternative_folder_name']
	DBNAME: str = CONFIG['database_name']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	top_dir = f'{WORKING_DIR}/{ALT_NAME.replace(" ", "_")}/repeatmasking' if ALT_NAME else f'{WORKING_DIR}/{SPECIES_NAME.replace(" ", "_")}/repeatmasking'
	top_out = f'{OUTPUT_DIR}/{ALT_NAME.replace(" ", "_")}/repeatmasking' if ALT_NAME else f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/repeatmasking'

	build_database = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_build_database',
		template=build_repeatmodeller_database(
			genome_assembly_file=GENOME_ASSEMBLY,
			output_directory=top_dir,
			species_name=SPECIES_NAME
		)
	)
	
	repmod = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_RepeatModeler',
		template=repeatmodeler(
			database=build_database.outputs['db_files'],
			output_directory=top_dir,
			species_name=SPECIES_NAME
		)
	)

	repeatmasker_database = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_RepeatMasker_{DBNAME}',
		template=repeatmasker(
			genome_assembly_file=GENOME_ASSEMBLY,
			library_file=REPEAT_DATABASE,
			output_directory=top_dir,
			run_name=f'{DBNAME}_run'
		)
	)

	repeatmasker_repmod = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_RepeatMasker_repmod',
		template=repeatmasker(
			genome_assembly_file=repeatmasker_database.outputs['repmaskout'][1],
			library_file=repmod.outputs['all'],
			output_directory=top_dir,
			run_name=f'RepMod_run'
		)
	)

	combine = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_combine_RepeatMasker_output',
		template=combine_repeatmasker_runs(
			repeatmasker_run1=repeatmasker_database.outputs['repmaskout'],
			repeatmasker_run2=repeatmasker_repmod.outputs['repmaskout'],
			library_file=REPEAT_DATABASE,
			output_directory=top_dir,
			species_name=SPECIES_NAME
		)
	)

	softmasking = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_softmasking',
		template=mask_assembly(
			genome_assembly_file=GENOME_ASSEMBLY,
			annotation_file=combine.outputs['gff'],
			output_directory=top_out if OUTPUT_DIR else top_dir,
			species_name=SPECIES_NAME
		)
	)

	return gwf