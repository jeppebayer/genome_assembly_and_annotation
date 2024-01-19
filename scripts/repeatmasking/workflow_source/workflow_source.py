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
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    GENOME_ASSEMBLY: str = config['genome_assembly']
    REPEAT_DATABASE: str = config['repeat_database']
    WORKING_DIR: str = config['working_directory_path']
    OUTPUT_DIR: str = config['output_directory_path']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{WORKING_DIR}/{SPECIES_NAME.replace(" ", "_")}/repeatmasking'

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
            working_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    repeatmasker_database = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_RepeatMasker_repbase_arthropoda',
        template=repeatmasker(
            genome_assembly_file=GENOME_ASSEMBLY,
            library_file=REPEAT_DATABASE,
            output_directory=top_dir,
            run_name='RepBase_arthropoda_run'
        )
    )

    repeatmasker_repmod = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_RepeatMasker_repmod',
        template=repeatmasker(
            genome_assembly_file=repeatmasker_database.outputs['repmaskout'][1],
            library_file=repmod.outputs['all'],
            output_directory=top_dir,
            run_name='RepMod_run'
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
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    return gwf