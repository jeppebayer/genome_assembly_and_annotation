#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def hic_scaffolding_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Scaffolds draft genome assembly using Hi-C data
    
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
    HIC_READS: list = config['hic_sequence_files']
    DRAFT_GENOME: str = config['draft_genome_file']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}'
    os.makedirs(top_dir, exist_ok=True)

    index = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_index',
        template=index_reference(
            reference_genome_file=DRAFT_GENOME
        )
    )
    
    hic_align = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_HiC_alignment',
        template=hic_alignment_to_draft_assembly(
            hic_seqeuence_files=HIC_READS,
            draft_genome_file=DRAFT_GENOME,
            reference_indices=index.outputs['bwa'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    mark_duplicates = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_mark_duplicates',
        template=mark_duplicates_picard(
            alignment_bam_file=hic_align.outputs['bam'],
            output_directory=top_dir
        )
    )

    hic_scaffolding = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_HiC_scaffolding',
        template=hic_scaffolding_yahs(
            draft_genome_file=DRAFT_GENOME,
            hic_to_draft_bam_file=mark_duplicates.outputs['markdup'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    conversion = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_juicer_conversion',
        template=yahs_conversion_manual_curation(
            hic_bin_file=hic_scaffolding.outputs['bin'],
            scaffolds_final_agp_file=hic_scaffolding.outputs['final_agp'],
            draft_assembly_fai_index_file=index.outputs['fai'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    matrix = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_matrix',
        template=contact_matrix_manual_curation(
            JBAT_text_file=conversion.outputs['jbat'][0],
            JBAT_log_file=conversion.outputs['jbat'][4],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    return gwf

def juicer_hic_scaffolding_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Scaffolds draft genome assembly using Hi-C data with the :script:`juicer` pipeline.
    
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
    HIC_READS: list = config['hic_sequence_files']
    DRAFT_GENOME: str = config['draft_genome_file']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}'
    os.makedirs(top_dir, exist_ok=True)

    setup = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_juicer_setup',
        template=setup_for_juicer(
            hic_sequence_files=HIC_READS,
            draft_genome_assembly_file=DRAFT_GENOME,
            output_directory=top_dir
        )
    )
    
    juicer = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_juicer_pipeline_start',
        template=juicer_pipeline(
            draft_genome_assembly_file=setup.outputs['reference'],
            chrom_sizes_file=setup.outputs['sizes'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    return gwf