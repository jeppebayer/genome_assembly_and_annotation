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
    DATABASE: int = config['database']
    PROTEIN_DB: str = config['protein_database_file']
    REFERENCE: str = config['referene_genome_file']
    GTF: str = config['gtf_genome_annotation_file']
    
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
    
    rna_alignment = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_STAR_alignment',
        template=star_alignment(
            rna_sequence_files=RNA_READS,
            star_index_directory=f'{top_dir}/indices',
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    if DATABASE == 0:
        braker_no_db

    elif DATABASE == 1:
        braker_db_exists = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_braker3',
            template=braker3(
                genome_assembly_file=GENOME_ASSEMBLY,
                rna_alignment_bam=rna_alignment.outputs['bam'],
                protein_database_file=PROTEIN_DB,
                output_directory=top_dir,
                species_name=SPECIES_NAME
            )
        )
    elif DATABASE == 2:
        make_database = gwf.target_from_template(
            name=f'Create_protein_database',
            template=make_protein_db(
                reference_genome_file=REFERENCE,
                gtf_annotation_file=GTF,
                output_directory=OUTPUT_DIR
            )
        )
        braker_db_exists_not = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_braker3',
            template=braker3(
                genome_assembly_file=GENOME_ASSEMBLY,
                rna_alignment_bam=rna_alignment.outputs['bam'],
                protein_database_file=make_database.outputs['proteins'],
                output_directory=top_dir,
                species_name=SPECIES_NAME
            )
        )

    return gwf