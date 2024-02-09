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

#------------------------------------------------------------------------
########################## Draft assembly HiFi ##########################
#------------------------------------------------------------------------

def hifiadapterfilt(pacbio_hifi_file: str, output_directory: str, hifiadapterfilt_directory: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/HiFiAdapterFilt'):
	"""
	Template: Removes remaining adapters from PacBio HiFi reads.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	os.makedirs(f'{output_directory}/draft_assembly/HiFiAdapterFilt', exist_ok=True)
	working_dir = f'{output_directory}/draft_assembly/HiFiAdapterFilt'
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
	
	cp {pacbio_hifi_file} {output_directory}/draft_assembly/HiFiAdapterFilt/prog{ext}

	export PATH=$PATH:{hifiadapterfilt_directory}
	export PATH=$PATH:{hifiadapterfilt_directory}/DB

	bash {hifiadapterfilt_directory}/pbadapterfilt.sh \
		-t {options['cores']} \
		-p prog \
		-o {output_directory}/draft_assembly/HiFiAdapterFilt

	mv {output_directory}/draft_assembly/HiFiAdapterFilt/prog.filt.fastq.gz {outputs['filt']}
    mv {output_directory}/draft_assembly/HiFiAdapterFilt/prog.contaminant.blastout {outputs['cont']}
    mv {output_directory}/draft_assembly/HiFiAdapterFilt/prog.blocklist {outputs['block']}
    mv {output_directory}/draft_assembly/HiFiAdapterFilt/prog.stats {outputs['stats']}
    rm {output_directory}/draft_assembly/HiFiAdapterFilt/prog{ext}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def hifiasm_primary(hifi_sequence_file: str, output_directory: str, species_name: str, similarity_threshold: int = 0.1, purge_level: int = 3):
	"""
	Template: Create draft genome assembly using PacBio HiFi data.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hifi': hifi_sequence_file}
	outputs = {'fasta': f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_ctg.fasta',
			   'primary': [f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_ctg.gfa',
				  		   f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_ctg.lowQ.bed',
						   f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_ctg.noseq.gfa'],
			   'alt': [f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.a_ctg.gfa',
				  	   f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.a_ctg.lowQ.bed',
					   f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.a_ctg.noseq.gfa'],
			   'raw': [f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.r_utg.gfa',
				  	   f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.r_utg.lowQ.bed',
					   f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.r_utg.noseq.gfa'],
			   'nobub': [f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_utg.gfa',
				  	   	 f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_utg.lowQ.bed',
					   	 f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_utg.noseq.gfa'],
			   'other': [f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.ec.bin',
			   			 f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.ovlp.reverse.bin',
						 f'{output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.ovlp.source.bin']}
	protect = [outputs['fasta'], outputs['primary'][0], outputs['primary'][1], outputs['primary'][2]]
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
	
	[ -d {output_directory}/draft_assembly/hifiasm_primary ] || mkdir -p {output_directory}/draft_assembly/hifiasm_primary
	
	hifiasm \
		-t {options['cores']} \
		-s {similarity_threshold} \
		-l {purge_level} \
		--primary \
		-o {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog \
		{hifi_sequence_file}
	
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.p_ctg.gfa {outputs['primary'][0]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.p_ctg.lowQ.bed {outputs['primary'][1]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.p_ctg.noseq.gfa {outputs['primary'][2]}
    mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.a_ctg.gfa {outputs['alt'][0]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.a_ctg.lowQ.bed {outputs['alt'][1]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.a_ctg.noseq.gfa {outputs['alt'][2]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.r_utg.gfa {outputs['raw'][0]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.r_utg.lowQ.bed {outputs['raw'][1]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.r_utg.noseq.gfa {outputs['raw'][2]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.p_utg.gfa {outputs['nobub'][0]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.p_utg.lowQ.bed {outputs['nobub'][1]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.p_utg.noseq.gfa {outputs['nobub'][2]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.ec.bin {outputs['other'][0]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.ovlp.reverse.bin {outputs['other'][1]}
	mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.prog.ovlp.source.bin {outputs['other'][2]}
	
	awk \
        'BEGIN{{FS="\\t"}}
        {{if ($0 ~ /^S/)
            {{print ">"$2"\\n"$3}}
		}}' \
        {outputs['primary'][0]} \
    | fold \
        > {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_ctg.prog.fasta

    mv {output_directory}/draft_assembly/hifiasm_primary/{species_abbreviation(species_name)}.p_ctg.prog.fasta {outputs['fasta']}

	awk \
		'BEGIN{{OFS="\\t"}}
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
		> {output_directory}/draft_assembly/hifiasm_primary/sequences.tsv

	last_line=($(tail -n 1 {output_directory}/draft_assembly/hifiasm_primary/sequences.tsv))

	mv {output_directory}/draft_assembly/hifiasm_primary/sequences.tsv {output_directory}/draft_assembly/hifiasm_primary/sequences_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def hifiasm_hic(hifi_sequence_file: str, hic_sequence_files: list, output_directory: str, species_name: str, similarity_threshold: int = 0.1, purge_level: int = 3):
	"""
	Template: Create draft genome assembly using PacBio HiFi data.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hifi': hifi_sequence_file,
			  'hic': hic_sequence_files}
	outputs = {'fasta': [f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_ctg.fasta',
					  	 f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap1.p_ctg.fasta',
						 f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap2.p_ctg.fasta'],
			   'primary': [f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_ctg.gfa',
				  		   f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_ctg.lowQ.bed',
						   f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_ctg.noseq.gfa'],
			   'hap1': [f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap1.p_ctg.gfa',
				  		f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap1.p_ctg.lowQ.bed',
						f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap1.p_ctg.noseq.gfa'],
			   'hap2': [f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap2.p_ctg.gfa',
				  		f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap2.p_ctg.lowQ.bed',
						f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap2.p_ctg.noseq.gfa'],
			   'raw': [f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.r_utg.gfa',
				  	   f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.r_utg.lowQ.bed',
					   f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.r_utg.noseq.gfa'],
			   'nobub': [f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_utg.gfa',
				  	   	 f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_utg.lowQ.bed',
					   	 f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_utg.noseq.gfa'],
			   'other': [f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.ec.bin',
			   			 f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.ovlp.reverse.bin',
						 f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.ovlp.source.bin',
						 f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.lk.bin',
						 f'{output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.tlb.bin']}
	protect = outputs['fasta'] + outputs['hap1'] + outputs['hap2']
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
	
	[ -d {output_directory}/draft_assembly/hifiasm_hic ] || mkdir -p {output_directory}/draft_assembly/hifiasm_hic
	
	hifiasm \
		-t {options['cores']} \
		-s {similarity_threshold} \
		-l {purge_level} \
		-o {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog \
		--h1 {hic_sequence_files[0]} \
		--h2 {hic_sequence_files[1]} \
		{hifi_sequence_file}
	
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.p_ctg.gfa {outputs['primary'][0]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.p_ctg.lowQ.bed {outputs['primary'][1]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.p_ctg.noseq.gfa {outputs['primary'][2]}
    mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.hap1.p_ctg.gfa {outputs['hap1'][0]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.hap1.p_ctg.lowQ.bed {outputs['hap1'][1]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.hap1.p_ctg.noseq.gfa {outputs['hap1'][2]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.hap2.p_ctg.gfa {outputs['hap2'][0]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.hap2.p_ctg.lowQ.bed {outputs['hap2'][1]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.hap2.p_ctg.noseq.gfa {outputs['hap2'][2]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.r_utg.gfa {outputs['raw'][0]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.r_utg.lowQ.bed {outputs['raw'][1]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.r_utg.noseq.gfa {outputs['raw'][2]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.p_utg.gfa {outputs['nobub'][0]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.p_utg.lowQ.bed {outputs['nobub'][1]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.p_utg.noseq.gfa {outputs['nobub'][2]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.ec.bin {outputs['other'][0]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.ovlp.reverse.bin {outputs['other'][1]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.ovlp.source.bin {outputs['other'][2]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.lk.bin {outputs['other'][3]}
	mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.prog.hic.tlb.bin {outputs['other'][4]}
	
	awk \
        'BEGIN{{FS="\\t"}}
        {{if ($0 ~ /^S/)
            {{print ">"$2"\\n"$3}}}}' \
        {outputs['primary'][0]} \
    | fold \
        > {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_ctg.prog.fasta

    mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.p_ctg.prog.fasta {outputs['fasta'][0]}

	awk \
        'BEGIN{{FS="\\t"}}
        {{if ($0 ~ /^S/)
            {{print ">"$2"\\n"$3}}}}' \
        {outputs['hap1'][0]} \
    | fold \
        > {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap1.p_ctg.prog.fasta

    mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap1.p_ctg.prog.fasta {outputs['fasta'][1]}

	awk \
        'BEGIN{{FS="\\t"}}
        {{if ($0 ~ /^S/)
            {{print ">"$2"\\n"$3}}}}' \
        {outputs['hap2'][0]} \
    | fold \
        > {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap2.p_ctg.prog.fasta

    mv {output_directory}/draft_assembly/hifiasm_hic/{species_abbreviation(species_name)}.hic.hap2.p_ctg.prog.fasta {outputs['fasta'][2]}
	
	awk \
		'BEGIN{{OFS="\\t"}}
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
		{outputs['fasta'][0]} \
		> {output_directory}/draft_assembly/hifiasm_hic/sequences_primary.tsv

	last_line=($(tail -n 1 {output_directory}/draft_assembly/hifiasm_hic/sequences_primary.tsv))

	mv {output_directory}/draft_assembly/hifiasm_hic/sequences_primary.tsv {output_directory}/draft_assembly/hifiasm_hic/sequences_primary_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	awk \
		'BEGIN{{OFS="\\t"}}
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
		{outputs['fasta'][1]} \
		> {output_directory}/draft_assembly/hifiasm_hic/sequences_hap1.tsv

	last_line=($(tail -n 1 {output_directory}/draft_assembly/hifiasm_hic/sequences_hap1.tsv))

	mv {output_directory}/draft_assembly/hifiasm_hic/sequences_hap1.tsv {output_directory}/draft_assembly/hifiasm_hic/sequences_hap1_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	awk \
		'BEGIN{{OFS="\\t"}}
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
		{outputs['fasta'][2]} \
		> {output_directory}/draft_assembly/hifiasm_hic/sequences_hap2.tsv

	last_line=($(tail -n 1 {output_directory}/draft_assembly/hifiasm_hic/sequences_hap2.tsv))

	mv {output_directory}/draft_assembly/hifiasm_hic/sequences_hap2.tsv {output_directory}/draft_assembly/hifiasm_hic/sequences_hap2_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def busco_genome(genome_assembly_file: str, busco_dataset: str, busco_download_path: str = '/faststorage/project/EcoGenetics/BACKUP/database/busco'):
	"""
	Template: Runs BUSCO analysis on genome assembly.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': genome_assembly_file}
	outputs = {'stattxt': f'{os.path.dirname(genome_assembly_file)}/busco_{os.path.basename(genome_assembly_file)}/short_summary.specific.{busco_dataset}.busco_{os.path.basename(genome_assembly_file)}.txt',
			   'statjson': f'{os.path.dirname(genome_assembly_file)}/busco_{os.path.basename(genome_assembly_file)}/short_summary.specific.{busco_dataset}.busco_{os.path.basename(genome_assembly_file)}.json'}
	protect = [outputs['stattxt'], outputs['statjson']]
	options = {
		'cores': 30,
		'memory': '50g',
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
	
	busco \
		--cpu {options['cores']} \
		--force \
		--in {genome_assembly_file} \
		--mode genome \
		--out busco_{os.path.basename(genome_assembly_file)} \
		--out_path {os.path.dirname(genome_assembly_file)} \
		--download_path {busco_download_path} \
		--lineage {busco_download_path}/lineages/{busco_dataset} \
		--tar \
		--offline
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def busco_protein(protein_sequence_file: str, busco_dataset: str, busco_download_path: str = '/faststorage/project/EcoGenetics/BACKUP/database/busco'):
	"""
	Template: Runs BUSCO analysis on protein sequences from an annotated gene set.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': protein_sequence_file}
	outputs = {'stattxt': f'{os.path.dirname(protein_sequence_file)}/busco_{os.path.basename(protein_sequence_file)}/short_summary.specific.{busco_dataset}.busco_{os.path.basename(protein_sequence_file)}.txt',
			   'statjson': f'{os.path.dirname(protein_sequence_file)}/busco_{os.path.basename(protein_sequence_file)}/short_summary.specific.{busco_dataset}.busco_{os.path.basename(protein_sequence_file)}.json'}
	protect = [outputs['stattxt'], outputs['statjson']]
	options = {
		'cores': 30,
		'memory': '50g',
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
	
	busco \
		--cpu {options['cores']} \
		--force \
		--in {protein_sequence_file} \
		--mode proteins \
		--out busco_{os.path.basename(protein_sequence_file)} \
		--out_path {os.path.dirname(protein_sequence_file)} \
		--download_path {busco_download_path} \
		--lineage {busco_download_path}/lineages/{busco_dataset} \
		--tar
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

########################## Purge duplicates ##########################
def purge_dups_1_map_hifi_to_genome(gemone_assembly_file: str, hifi_sequence_file: str, output_directory: str, species_name: str, directory_addition: str = None, round_number: int = 1):
	"""
	Template: Align PacBio HiFi sequence data to genome.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	if directory_addition:
		directory_addition = '_' + directory_addition
	else:
		directory_addition = ''
	inputs = {'genome': gemone_assembly_file,
		   	  'hifi': hifi_sequence_file}
	outputs = {'paf': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.paf.gz',
			   'stat': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/PB.stat',
			   'cov': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/PB.base.cov'}
	options = {
		'cores': 30,
		'memory': '60g',
		'walltime': '04:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments ] || mkdir -p {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments
	
	minimap2 \
		-x map-hifi \
		-t {options['cores']} \
		{gemone_assembly_file} \
		{hifi_sequence_file} \
	| gzip \
		-c \
		- \
		> {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.prog.paf.gz
	
	pbcstat \
		-O {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/ \
		{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.prog.paf.gz
	
	mv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.prog.paf.gz {outputs['paf']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def purge_dups_2_map_to_self(genome_assembly_file: str, output_directory: str, species_name: str, directory_addition: str = None, round_number: int = 1):
	"""
	Template: Splits assembly into equal pieces and maps it to itself.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	if directory_addition:
		directory_addition = '_' + directory_addition
	else:
		directory_addition = ''
	inputs = {'genome': genome_assembly_file}
	outputs = {'split': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.split.fasta',
			   'paf': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.self.paf.gz'}
	options = {
		'cores': 30,
		'memory': '60g',
		'walltime': '04:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments ] || mkdir -p {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments
	
	split_fa \
		{genome_assembly_file} \
		> {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.split.prog.fasta
	
	mv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.split.prog.fasta {outputs['split']}

	minimap2 \
		-x asm5 \
		-t {options['cores']} \
		-D \
		-P \
		{outputs['split']} \
		{outputs['split']} \
	| gzip \
		-c \
		- \
		> {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.self.prog.paf.gz
	
	mv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/alignments/{species_abbreviation(species_name)}.self.prog.paf.gz {outputs['paf']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def purge_dups_3_purge_duplicates(pb_stat_file: str, pb_base_cov_file: str, self_alignment_paf: str, genome_assembly_file: str, output_directory: str, species_name: str, directory_addition: str = None, round_number: int = 1):
	"""
	Template: Purge haplotigs and opverlaps using :script:`purge_dups`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	if directory_addition:
		directory_addition = '_' + directory_addition
	else:
		directory_addition = ''
	inputs = {'stat': pb_stat_file,
		   	  'cov': pb_base_cov_file,
			  'self': self_alignment_paf,
			  'genome': genome_assembly_file}
	outputs = {'cutoffs': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/cutoffs',
			   'dups': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/dups.bed',
			   'purged': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/{species_abbreviation(species_name)}.purged.fa',
			   'hap': f'{output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/{species_abbreviation(species_name)}.hap.fa'}
	options = {
		'cores': 2,
		'memory': '20g',
		'walltime': '04:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02} ] || mkdir -p {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}
	
	calcuts \
		{pb_stat_file} \
		> {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/cutoffs.prog
	
	mv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/cutoffs.prog {outputs['cutoffs']}

	purge_dups \
		-2 \
		-T {outputs['cutoffs']} \
		-c {pb_base_cov_file} \
		{self_alignment_paf} \
		> {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/dups.prog.bed

	mv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/dups.prog.bed {outputs['dups']}

	get_seqs \
		-p {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/{species_abbreviation(species_name)}.prog \
		{outputs['dups']} \
		{genome_assembly_file}

	mv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/{species_abbreviation(species_name)}.prog.purged.fa {outputs['purged']}
	mv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/{species_abbreviation(species_name)}.prog.hap.fa {outputs['hap']}
	
	awk \
		'BEGIN{{OFS="\\t"}}
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
		{outputs['purged']} \
		> {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/sequences.tsv

	last_line=($(tail -n 1 {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/sequences.tsv))

	mv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/sequences.tsv {output_directory}/draft_assembly/purge_dups{directory_addition}/{round_number:02}/sequences_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)