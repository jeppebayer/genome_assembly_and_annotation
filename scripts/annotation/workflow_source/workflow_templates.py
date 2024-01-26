#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

########################## Braker3 ##########################

def star_index(genome_assembly_file: str, output_directory: str):
	"""
	Template: Indexes genome assembly prior to :script:`STAR` alignment of RNA data.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	working_dir = output_directory
	inputs = {'assembly': genome_assembly_file}
	outputs = {}
	options = {
		'cores': 30,
		'memory': '240g',
		'walltime': '48:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate annotation
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/indices ] || mkdir -p {output_directory}/indices
	[ -d {output_directory}/tmp ] && rm -rf {output_directory}/tmp
	[ "$(ls -A {output_directory}/indices)" ] || rm -rf {output_directory}/indices/*
	
	STAR \
		--runThreadN {options['cores']} \
		--runMode genomeGenerate \
		--genomeDir {output_directory}/indices \
		--genomeFastaFiles {genome_assembly_file} \
		--outTmpDir {output_directory}/tmp
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)

def star_alignment(rna_sequence_files: list, star_index_directory: str, output_directory: str, species_name: str):
	"""
	Template: Align RNA data to indexed genome assembly using :script:`STAR`
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	working_dir = output_directory
	inputs = {}
	outputs = {}
	options = {
		'cores': 30,
		'memory': '240g',
		'walltime': '48:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate annotation
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/rna_alignment ] || mkdir -p {output_directory}/rna_alignment
	[ -d {output_directory}/tmp ] && rm -rf {output_directory}/tmp
	
	STAR \
		--runThreadsN {options['cores']} \
		--runMode alignReads \
		--genomeDir {star_index_directory} \
		--readFilesIn {" ".join(rna_sequence_files)} \
		--readFilesCommand zcat \
		--outFileNamePrefix {output_directory}/rna_alignment/{species_abbreviation(species_name)} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMstrandfield intronMotif \
		--outSAMattrRGline ID:{species_abbreviation(species_name)}_RNA SM:{species_abbreviation(species_name)}_RNA \
		--outTempDir {output_directory}/tmp
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)