#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
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
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    OUTPUT_DIR: str = config['output_directory_path']
    HIFI_SEQ: str = config['pacbio_hifi_sequence_file']
    HIC_READS: list = config['hic_sequence_files']
    BUSCO_DATASET: str = config['busco_dataset_name']
    PURGE_ROUNDS: int = config['purge_dups_rounds']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}'
    os.makedirs(top_dir, exist_ok=True)

    remove_adapters = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_adapterremoval',
        template=hifiadapterfilt(
            pacbio_hifi_file=HIFI_SEQ,
            output_directory=top_dir,
        )
    )
    
    # if HIC_READS:

    hifiasm_1 = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_hifiasm_primary',
        template=hifiasm_primary(
            hifi_sequence_file=remove_adapters.outputs['filt'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    hifiasm_2 = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_hifiasm_hic',
        template=hifiasm_hic(
            hifi_sequence_file=HIFI_SEQ,
            hic_sequence_files=HIC_READS,
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    busco_1 = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_busco_primary',
        template=busco_genome(
            genome_assembly_file=hifiasm_1.outputs['fasta'],
            busco_dataset=BUSCO_DATASET
        )
    )

    busco_2 = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_busco_hic',
        template=busco_genome(
            genome_assembly_file=hifiasm_2.outputs['fasta'][0],
            busco_dataset=BUSCO_DATASET
        )
    )

    for i in range(1, PURGE_ROUNDS + 1):
        if i == 1:
            in_file = hifiasm_2.outputs['fasta'][0]
        else:
            in_file = step3.outputs['purged']
        step1 = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_purge_dups_step_1_round_{i:02}',
            template=purge_dups_1_map_hifi_to_genome(
                gemone_assembly_file=in_file,
                hifi_sequence_file=HIFI_SEQ,
                output_directory=top_dir,
                species_name=SPECIES_NAME,
                round_number=i
            )
        )
        step2 = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_purge_dups_step_2_round_{i:02}',
            template=purge_dups_2_map_to_self(
                genome_assembly_file=in_file,
                output_directory=top_dir,
                species_name=SPECIES_NAME,
                round_number=i
            )
        )
        step3 = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_purge_dups_step_3_round_{i:02}',
            template=purge_dups_3_purge_duplicates(
                pb_stat_file=step1.outputs['stat'],
                pb_base_cov_file=step1.outputs['cov'],
                self_alignment_paf=step2.outputs['paf'],
                genome_assembly_file=in_file,
                output_directory=top_dir,
                species_name=SPECIES_NAME,
                round_number=i
            )
        )
        purge_busco = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_purge_dups_busco_round_{i:02}',
            template=busco_genome(
                genome_assembly_file=step3.outputs['purged'],
                busco_dataset=BUSCO_DATASET
            )
        )

    return gwf