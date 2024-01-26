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
	[ "$(ls -A {output_directory}/indices)" ] || rm -rf {output_directory}/indices/*
	
	STAR \
		--runThreadN {options['cores']} \
		--runMode genomeGenerate \
		--genomeDir {output_directory}/indices \
		--genomeFastaFiles {genome_assembly_file} \
		--sjbOverhang 149
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def star_alignment(rna_sequence_files: list, star_index_directory: str, output_directory: str, species_name: str):
	"""
	Template: Align RNA data to indexed genome assembly using :script:`STAR`
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
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
	
	STAR \
		--runThreadsN {options['cores']} \
		--genomeDir {star_index_directory} \
		--readFilesIn {" ".join(rna_sequence_files)} \
		--readFilesCommand zcat \
		--outFileNamePrefix {output_directory}/rna_alignment
	
	mv
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)