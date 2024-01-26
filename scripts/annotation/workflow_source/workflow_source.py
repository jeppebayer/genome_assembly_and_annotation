#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def braker3_annotation_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Annotates genome assembly using BRAKER3.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    OUTPUT_DIR: str = config['output_directory_path']
    GENOME_ASSEMBLY: str = config['genome_assembly_file']
    RNA_READS: list = config['rna_sequence_files']
    
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/annotation'
    os.makedirs(top_dir, exist_ok=True)

    index = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_STAR_index',
        template=star_index(
            genome_assembly_file=GENOME_ASSEMBLY,
            output_directory=top_dir
        )
    )
    
    return gwf