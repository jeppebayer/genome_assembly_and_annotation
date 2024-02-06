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

########################## Draft assembly HiFi ##########################

def hifiadapterfilt(pacbio_hifi_file: str, output_directory: str, hifiadapterfilt_directory: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/HiFiAdapterFilt'):
	"""
	Template: Removes remaining adapters from PacBio HiFi reads.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'pbhifi': pacbio_hifi_file}
	if pacbio_hifi_file.endswith('.gz'):
		prefix = os.path.splitext(os.path.splitext(os.path.basename(pacbio_hifi_file))[0])[0]
		ext = f'{os.path.splitext(os.path.splitext(os.path.basename(pacbio_hifi_file))[0])[1]}.gz'
	else:
		prefix = os.path.splitext(os.path.basename(pacbio_hifi_file))[0]
		ext = os.path.splitext(os.path.basename(pacbio_hifi_file))[1]
	outputs = {'filt': f'{output_directory}/draft_assembly/HiFiAdapterFilt/{prefix}.filt.fastq.gz',
               'cont': f'{output_directory}/draft_assembly/HiFiAdapterFilt/{prefix}.contaminant.blastout',
               'block': f'{output_directory}/draft_assembly/HiFiAdapterFilt/{prefix}.blocklist',
               'stats': f'{output_directory}/draft_assembly/HiFiAdapterFilt/{prefix}.stats'}
	protect = [outputs['filt'], outputs['stats']]
	options = {
		'cores': 30,
		'memory': '30g',
		'walltime': '02:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/draft_assembly/HiFiAdapterFilt ] || mkdir -p {output_directory}/draft_assembly/HiFiAdapterFilt
	
	export PATH=$PATH:{hifiadapterfilt_directory}
	export PATH=$PATH:{hifiadapterfilt_directory}/DB
	
	ln -s {pacbio_hifi_file} {output_directory}/prog.{ext}

	bash {hifiadapterfilt_directory}/pbadapterfilt.sh \
		-t {options['cores']} \
		-p {output_directory}/draft_assembly/HiFiAdapterFilt/prog \
		-o {output_directory}/draft_assembly/HiFiAdapterFilt

	mv {output_directory}HiFiAdapterFilt/prog.filt.fastq.gz {outputs['filt']}
    mv {output_directory}HiFiAdapterFilt/prog.contaminant.blastout {outputs['cont']}
    mv {output_directory}HiFiAdapterFilt/prog.blocklist {outputs['block']}
    mv {output_directory}HiFiAdapterFilt/prog.stats {outputs['stats']}
    rm {output_directory}/draft_assembly/HiFiAdapterFilt//prog.{ext}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def hifiasm_primary(hifi_sequence_file: str, output_directory: str, species_name: str, similarity_threshold: int = 0.1, purge_level: int = 3):
	"""
	Template: Create draft genome assembly using PacBio HiFi data
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hifi': hifi_sequence_file}
	outputs = {'fasta'}
	options = {
		'cores': 32,
		'memory': '100g',
		'walltime': '10:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/draft_assembly/hifiasm ] || mkdir -p {output_directory}/draft_assembly/hifiasm
	
	hifiasm \
		-t {options['cores']} \
		-s {similarity_threshold} \
		-l {purge_level} \
		--primary \
		-o {output_directory}/draft_assembly/hifiasm/{species_abbreviation(species_name)}.prog
	
	mv {output_directory}/draft_assembly/hifiasm/{species_abbreviation(species_name)}.asm.bp.p_ctg.gfa {outputs['primary']}
    mv {output_directory}/draft_assembly/hifiasm/{species_abbreviation(species_name)}.asm.a_ctg.gfa
    mv {output_directory}/draft_assembly/hifiasm/{species_abbreviation(species_name)}.asm.r_utg.gfa
    mv {output_directory}/draft_assembly/hifiasm/{species_abbreviation(species_name)}.asm.p_ctg.gfa
	
	awk \
        'BEGIN{{FS="\t"}}
        {{if ($0 ~ /^S/)
            {{print ">"$2"\n"$3}}}}' \
        {outputs['primary']} \
    | fold \
        > {output_directory}/draft_assembly/hifiasm/{species_abbreviation(species_name)}.bp.p_ctg.prog.fasta

    mv {output_directory}/draft_assembly/hifiasm/{species_abbreviation(species_name)}.bp.p_ctg.prog.fasta {outputs['fasta']}

	awk \
		'BEGIN{{OFS="\t"}}
		{{if ($0 ~ /^>/)
			{{if (sequence_length)
				{{print sequence_name, sequence_length}}
			total_length += sequence_length
			sequence_name = $0
			sequence_length = 0
			sequence_number += 1
			next
			}}
		sequence_length += length($0)
		}}
		END{{if (sequence_length)
			{{print sequence_name, sequence_length}}
		print sequence_number"_sequences", total_length + sequence_length
		}}' \
		{outputs['fasta']} \
		> {output_directory}/draft_assembly/hifiasm/sequences.tsv

	last_line=($(tail -n 1 {output_directory}/draft_assembly/hifiasm/sequences.tsv))

	mv {output_directory}/draft_assembly/hifiasm/sequences.tsv {output_directory}/draft_assembly/hifiasm/sequences_"${{last_line[0]%_*}}"_"${{last_line[1]}}.tsv

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)