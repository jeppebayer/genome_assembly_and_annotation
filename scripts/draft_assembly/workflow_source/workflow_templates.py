#!/bin/env python3
from gwf import AnonymousTarget
from gwf.executors import Conda
import os, yaml

########################## Functions ##########################

def speciesAbbreviation(species_name: str) -> str:
	"""Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

def software_versions(environments: list, packages: list) -> dict:
	environmentDict = {}
	for environment in environments:
		environmentDict[environment] = []
		if os.path.isdir(environment):
			environmentJson = yaml.safe_load(os.popen(f'conda list --prefix {environment} --json'))
			for package in environmentJson:
				if package['name'] in packages:
					environmentDict[environment].append({'name': package['name'], 'version': package['version']})
		elif type(environment) == str:
			environmentJson = yaml.safe_load(os.popen(f'conda list --name {environment} --json'))
			for package in environmentJson:
				if package['name'] in packages:
					environmentDict[environment].append({'name': package['name'], 'version': package['version']})
	return environmentDict

def software_versions_to_string(environmentDict: dict) -> str:
	outputString = '#environment\tpackage\tversion\n'
	for environment in environmentDict:
		for package in environmentDict[environment]:
			outputString += f"{environment}\t{package['name']}\t{package['version']}\n"
	return outputString

#------------------------------------------------------------------------
########################## Draft assembly HiFi ##########################
#------------------------------------------------------------------------

# def hifiadapterfilt(pacbioHifiFiles: list, outputDirectory: str, speciesName: str, hifiadapterfiltDirectory: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/HiFiAdapterFilt'):
# 	"""
# 	Template: Removes remaining adapters from PacBio HiFi reads.
	
# 	Template I/O::
	
# 		inputs = {}
# 		outputs = {}
	
# 	:param
# 	"""
# 	inputs = {'pbhifi': pacbioHifiFiles}
# 	if pacbioHifiFiles[0].endswith('.gz'):
# 		infile = f'<(zcat {" ".join(pacbioHifiFiles)})'
# 	else:
# 		infile = f'<(cat {" ".join(pacbioHifiFiles)})'
# 	outputs = {'filt': f'{outputDirectory}/HiFiAdapterFilt/{speciesAbbreviation(speciesName)}.filt.fastq.gz',
#                'cont': f'{outputDirectory}/HiFiAdapterFilt/{speciesAbbreviation(speciesName)}.contaminant.blastout',
#                'block': f'{outputDirectory}/HiFiAdapterFilt/{speciesAbbreviation(speciesName)}.blocklist',
#                'stats': f'{outputDirectory}/HiFiAdapterFilt/{speciesAbbreviation(speciesName)}.stats'}
# 	protect = [outputs['filt'], outputs['stats']]
# 	options = {
# 		'cores': 30,
# 		'memory': '30g',
# 		'walltime': '06:00:00'
# 	}
# 	spec = f"""
# 	# Sources environment
# 	if [ "$USER" == "jepe" ]; then
# 		source /home/"$USER"/.bashrc
# 		source activate assembly
# 	fi
	
# 	echo "START: $(date)"
# 	echo "JobID: $SLURM_JOBID"
	
# 	[ -d {outputDirectory}/HiFiAdapterFilt ] || mkdir -p {outputDirectory}/HiFiAdapterFilt

# 	cd {outputDirectory}/HiFiAdapterFilt

# 	bgzip \\
# 		-c \\
# 		{infile} \\
# 		> prog.fastq.gz

# 	export PATH=$PATH:{hifiadapterfiltDirectory}
# 	export PATH=$PATH:{hifiadapterfiltDirectory}/DB

# 	bash {hifiadapterfiltDirectory}/pbadapterfilt.sh \\
# 		-t {options['cores']} \\
# 		-p prog \\
# 		-o {outputDirectory}/HiFiAdapterFilt

# 	mv {outputDirectory}/HiFiAdapterFilt/prog.filt.fastq.gz {outputs['filt']}
#     mv {outputDirectory}/HiFiAdapterFilt/prog.contaminant.blastout {outputs['cont']}
#     mv {outputDirectory}/HiFiAdapterFilt/prog.blocklist {outputs['block']}
#     mv {outputDirectory}/HiFiAdapterFilt/prog.stats {outputs['stats']}
#     rm {outputDirectory}/HiFiAdapterFilt/prog.fastq.gz
	
# 	echo "END: $(date)"
# 	echo "$(jobinfo "$SLURM_JOBID")"
# 	"""
# 	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def hifiadapterfilt(pacbioHifiFiles: list, outputDirectory: str, speciesName: str, environment: str, hifiAdapterFiltDirectory: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/HiFiAdapterFilt'):
	"""
	Template: Removes remaining adapters from PacBio HiFi reads.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'pbhifi': pacbioHifiFiles}
	outputs = {'filt': f'{outputDirectory}/HiFiAdapterFilt/{speciesAbbreviation(speciesName)}.filt.fastq.gz',
               'cont': f'{outputDirectory}/HiFiAdapterFilt/{speciesAbbreviation(speciesName)}.contaminant.blastout',
               'block': f'{outputDirectory}/HiFiAdapterFilt/{speciesAbbreviation(speciesName)}.blocklist',
               'stats': f'{outputDirectory}/HiFiAdapterFilt/{speciesAbbreviation(speciesName)}.stats'}
	protect = [outputs['filt'], outputs['stats']]
	options = {
		'cores': 30,
		'memory': '30g',
		'walltime': '06:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/HiFiAdapterFilt ] || mkdir -p {outputDirectory}/HiFiAdapterFilt

	cd {outputDirectory}/HiFiAdapterFilt

	export PATH=$PATH:{hifiAdapterFiltDirectory}
	export PATH=$PATH:{hifiAdapterFiltDirectory}/DB

	bgzip \\
		-c \\
		{'<(zcat ' + ' '.join(pacbioHifiFiles) + ')' if pacbioHifiFiles[0].endswith('.gz') else '<(cat ' + ' '.join(pacbioHifiFiles) + ')'} \\
		> prog.fastq.gz

	bash {hifiAdapterFiltDirectory}/hifiadapterfilt.sh \\
		-t {options['cores']} \\
		-p prog \\
		-o {outputDirectory}/HiFiAdapterFilt

	mv {outputDirectory}/HiFiAdapterFilt/prog.filt.fastq.gz {outputs['filt']}
	mv {outputDirectory}/HiFiAdapterFilt/prog.contaminant.blastout {outputs['cont']}
	mv {outputDirectory}/HiFiAdapterFilt/prog.blocklist {outputs['block']}
	mv {outputDirectory}/HiFiAdapterFilt/prog.stats {outputs['stats']}
	rm {outputDirectory}/HiFiAdapterFilt/prog.fastq.gz
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec, executor=Conda(environment))

def hifiasm_primary(hifiSequenceFile: str, outputDirectory: str, speciesName: str, environment: str, similarityThreshold: float = 0.1, purgeLevel: int = 3):
	"""
	Template: Create draft genome assembly using PacBio HiFi data.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hifi': hifiSequenceFile}
	outputs = {'fasta': f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_ctg.fasta',
			   'primary': [f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_ctg.gfa',
				  		   f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_ctg.lowQ.bed',
						   f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_ctg.noseq.gfa'],
			   'alt': [f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.a_ctg.gfa',
				  	   f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.a_ctg.lowQ.bed',
					   f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.a_ctg.noseq.gfa'],
			   'raw': [f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.r_utg.gfa',
				  	   f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.r_utg.lowQ.bed',
					   f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.r_utg.noseq.gfa'],
			   'nobub': [f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_utg.gfa',
				  	   	 f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_utg.lowQ.bed',
					   	 f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_utg.noseq.gfa'],
			   'other': [f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.ec.bin',
			   			 f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.ovlp.reverse.bin',
						 f'{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.ovlp.source.bin']}
	protect = [outputs['fasta'], outputs['primary'][0], outputs['primary'][1], outputs['primary'][2]]
	options = {
		'cores': 32,
		'memory': '300g',
		'walltime': '72:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/hifiasm_primary ] || mkdir -p {outputDirectory}/hifiasm_primary
	
	hifiasm \\
		-t {options['cores']} \\
		-s {similarityThreshold} \\
		-l {purgeLevel} \\
		--primary \\
		-o {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog \\
		{hifiSequenceFile}
	
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.p_ctg.gfa {outputs['primary'][0]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.p_ctg.lowQ.bed {outputs['primary'][1]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.p_ctg.noseq.gfa {outputs['primary'][2]}
    mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.a_ctg.gfa {outputs['alt'][0]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.a_ctg.lowQ.bed {outputs['alt'][1]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.a_ctg.noseq.gfa {outputs['alt'][2]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.r_utg.gfa {outputs['raw'][0]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.r_utg.lowQ.bed {outputs['raw'][1]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.r_utg.noseq.gfa {outputs['raw'][2]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.p_utg.gfa {outputs['nobub'][0]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.p_utg.lowQ.bed {outputs['nobub'][1]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.p_utg.noseq.gfa {outputs['nobub'][2]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.ec.bin {outputs['other'][0]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.ovlp.reverse.bin {outputs['other'][1]}
	mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.prog.ovlp.source.bin {outputs['other'][2]}
	
	awk \\
        'BEGIN{{
			FS = "\\t"
		}}
        {{
			if ($0 ~ /^S/)
            {{
				print ">"$2"\\n"$3
			}}
		}}' \\
        {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_ctg.gfa \\
    | fold \\
        > {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_ctg.prog.fasta

    mv {outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_ctg.prog.fasta {outputs['fasta']}

	awk \\
		'BEGIN{{
			OFS="\\t"
		}}
		{{
			if ($0 ~ /^>/)
			{{
				if (sequence_length)
				{{
					print sequence_name, sequence_length
				}}
				total_length += sequence_length
				sequence_name = $0
				sequence_length = 0
				sequence_number += 1
				next
			}}
			sequence_length += length($0)
		}}
		END{{
			if (sequence_length)
			{{
				print sequence_name, sequence_length
			}}
			print sequence_number"_sequences", total_length + sequence_length
		}}' \\
		{outputDirectory}/hifiasm_primary/{speciesAbbreviation(speciesName)}.p_ctg.fasta \\
		> {outputDirectory}/hifiasm_primary/sequences.tsv

	last_line=($(tail -n 1 {outputDirectory}/hifiasm_primary/sequences.tsv))

	mv {outputDirectory}/hifiasm_primary/sequences.tsv {outputDirectory}/hifiasm_primary/sequences_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec, executor=Conda(environment))

def hifiasm_hic(hifiSequenceFile: str, hicSequenceFiles: list, outputDirectory: str, speciesName: str, environment: str, similarityThreshold: float = 0.1, purgeLevel: int = 3):
	"""
	Template: Create draft genome assembly using PacBio HiFi data.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hifi': hifiSequenceFile,
			  'hic': hicSequenceFiles}
	outputs = {'fasta': [f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_ctg.fasta',
					  	 f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap1.p_ctg.fasta',
						 f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap2.p_ctg.fasta'],
			   'primary': [f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_ctg.gfa',
				  		   f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_ctg.lowQ.bed',
						   f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_ctg.noseq.gfa'],
			   'hap1': [f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap1.p_ctg.gfa',
				  		f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap1.p_ctg.lowQ.bed',
						f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap1.p_ctg.noseq.gfa'],
			   'hap2': [f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap2.p_ctg.gfa',
				  		f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap2.p_ctg.lowQ.bed',
						f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap2.p_ctg.noseq.gfa'],
			   'raw': [f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.r_utg.gfa',
				  	   f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.r_utg.lowQ.bed',
					   f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.r_utg.noseq.gfa'],
			   'nobub': [f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_utg.gfa',
				  	   	 f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_utg.lowQ.bed',
					   	 f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_utg.noseq.gfa'],
			   'other': [f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.ec.bin',
			   			 f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.ovlp.reverse.bin',
						 f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.ovlp.source.bin',
						 f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.lk.bin',
						 f'{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.tlb.bin']}
	protect = outputs['fasta'] + outputs['hap1'] + outputs['hap2']
	options = {
		'cores': 32,
		'memory': '300g',
		'walltime': '72:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/hifiasm_hic ] || mkdir -p {outputDirectory}/hifiasm_hic
	
	hifiasm \\
		-t {options['cores']} \\
		-s {similarityThreshold} \\
		-l {purgeLevel} \\
		-o {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog \\
		--h1 {hicSequenceFiles[0]} \\
		--h2 {hicSequenceFiles[1]} \\
		{hifiSequenceFile}
	
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.p_ctg.gfa {outputs['primary'][0]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.p_ctg.lowQ.bed {outputs['primary'][1]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.p_ctg.noseq.gfa {outputs['primary'][2]}
    mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.hap1.p_ctg.gfa {outputs['hap1'][0]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.hap1.p_ctg.lowQ.bed {outputs['hap1'][1]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.hap1.p_ctg.noseq.gfa {outputs['hap1'][2]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.hap2.p_ctg.gfa {outputs['hap2'][0]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.hap2.p_ctg.lowQ.bed {outputs['hap2'][1]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.hap2.p_ctg.noseq.gfa {outputs['hap2'][2]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.r_utg.gfa {outputs['raw'][0]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.r_utg.lowQ.bed {outputs['raw'][1]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.r_utg.noseq.gfa {outputs['raw'][2]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.p_utg.gfa {outputs['nobub'][0]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.p_utg.lowQ.bed {outputs['nobub'][1]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.p_utg.noseq.gfa {outputs['nobub'][2]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.ec.bin {outputs['other'][0]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.ovlp.reverse.bin {outputs['other'][1]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.ovlp.source.bin {outputs['other'][2]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.lk.bin {outputs['other'][3]}
	mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.prog.hic.tlb.bin {outputs['other'][4]}
	
	awk \\
        'BEGIN{{
			FS = "\\t"
		}}
        {{
			if ($0 ~ /^S/)
            {{
				print ">"$2"\\n"$3
			}}
		}}' \\
        {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_ctg.gfa \\
    | fold \\
        > {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_ctg.prog.fasta

    mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_ctg.prog.fasta {outputs['fasta'][0]}

	awk \\
        'BEGIN{{
			FS = "\\t"
		}}
        {{
			if ($0 ~ /^S/)
            {{
				print ">"$2"\\n"$3
			}}
		}}' \\
        {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap1.p_ctg.gfa \\
    | fold \\
        > {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap1.p_ctg.prog.fasta

    mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap1.p_ctg.prog.fasta {outputs['fasta'][1]}

	awk \\
        'BEGIN{{
			FS = "\\t"
		}}
        {{
			if ($0 ~ /^S/)
            {{
				print ">"$2"\\n"$3
			}}
		}}' \\
        {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap2.p_ctg.gfa \\
    | fold \\
        > {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap2.p_ctg.prog.fasta

    mv {outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap2.p_ctg.prog.fasta {outputs['fasta'][2]}
	
	awk \\
		'BEGIN{{
			OFS = "\\t"
		}}
		{{
			if ($0 ~ /^>/)
			{{
				if (sequence_length)
				{{
					print sequence_name, sequence_length
				}}
				total_length += sequence_length
				sequence_name = $0
				sequence_length = 0
				sequence_number += 1
				next
			}}
			sequence_length += length($0)
		}}
		END{{
			if (sequence_length)
			{{
				print sequence_name, sequence_length
			}}
			print sequence_number"_sequences", total_length + sequence_length
		}}' \\
		{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.p_ctg.fasta \\
		> {outputDirectory}/hifiasm_hic/sequences_primary.tsv

	last_line=($(tail -n 1 {outputDirectory}/hifiasm_hic/sequences_primary.tsv))

	mv {outputDirectory}/hifiasm_hic/sequences_primary.tsv {outputDirectory}/hifiasm_hic/sequences_primary_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	awk \\
		'BEGIN{{
			OFS = "\\t"
		}}
		{{
			if ($0 ~ /^>/)
			{{
				if (sequence_length)
				{{
					print sequence_name, sequence_length
				}}
				total_length += sequence_length
				sequence_name = $0
				sequence_length = 0
				sequence_number += 1
				next
			}}
			sequence_length += length($0)
		}}
		END{{
			if (sequence_length)
			{{
				print sequence_name, sequence_length
			}}
			print sequence_number"_sequences", total_length + sequence_length
		}}' \\
		{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap1.p_ctg.fasta \\
		> {outputDirectory}/hifiasm_hic/sequences_hap1.tsv

	last_line=($(tail -n 1 {outputDirectory}/hifiasm_hic/sequences_hap1.tsv))

	mv {outputDirectory}/hifiasm_hic/sequences_hap1.tsv {outputDirectory}/hifiasm_hic/sequences_hap1_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	awk \\
		'BEGIN{{
			OFS = "\\t"
		}}
		{{
			if ($0 ~ /^>/)
			{{
				if (sequence_length)
				{{
					print sequence_name, sequence_length
				}}
				total_length += sequence_length
				sequence_name = $0
				sequence_length = 0
				sequence_number += 1
				next
			}}
			sequence_length += length($0)
		}}
		END{{
			if (sequence_length)
			{{
				print sequence_name, sequence_length
			}}
			print sequence_number"_sequences", total_length + sequence_length
		}}' \\
		{outputDirectory}/hifiasm_hic/{speciesAbbreviation(speciesName)}.hic.hap2.p_ctg.fasta \\
		> {outputDirectory}/hifiasm_hic/sequences_hap2.tsv

	last_line=($(tail -n 1 {outputDirectory}/hifiasm_hic/sequences_hap2.tsv))

	mv {outputDirectory}/hifiasm_hic/sequences_hap2.tsv {outputDirectory}/hifiasm_hic/sequences_hap2_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec, executor=Conda(environment))

########################## Quality Assessment ##########################
def busco_genome(genomeAssemblyFile: str, buscoDataset: str, environment: str, buscoDownloadPath: str = '/faststorage/project/EcoGenetics/databases/BUSCO'):
	"""
	Template: Runs BUSCO analysis on genome assembly.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': genomeAssemblyFile}
	outputs = {'stattxt': f'{os.path.dirname(genomeAssemblyFile)}/busco_{os.path.basename(genomeAssemblyFile)}/short_summary.specific.{buscoDataset}.busco_{os.path.basename(genomeAssemblyFile)}.txt',
			   'statjson': f'{os.path.dirname(genomeAssemblyFile)}/busco_{os.path.basename(genomeAssemblyFile)}/short_summary.specific.{buscoDataset}.busco_{os.path.basename(genomeAssemblyFile)}.json'}
	protect = [outputs['stattxt'], outputs['statjson']]
	options = {
		'cores': 30,
		'memory': '200g',
		'walltime': '10:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	cd {os.path.dirname(genomeAssemblyFile)}

	busco \\
		--cpu {options['cores']} \\
		--force \\
		--metaeuk \\
		--in {genomeAssemblyFile} \\
		--mode genome \\
		--out busco_{os.path.basename(genomeAssemblyFile)} \\
		--out_path {os.path.dirname(genomeAssemblyFile)} \\
		--download_path {buscoDownloadPath} \\
		--lineage {buscoDownloadPath}/lineages/{buscoDataset} \\
		--tar \\
		--offline

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec, executor=Conda(environment))

def busco_protein(proteinSequenceFile: str, buscoDataset: str, environment: str, buscoDownloadPath: str = '/faststorage/project/EcoGenetics/databases/BUSCO'):
	"""
	Template: Runs BUSCO analysis on protein sequences from an annotated gene set.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': proteinSequenceFile}
	outputs = {'stattxt': f'{os.path.dirname(proteinSequenceFile)}/busco_{os.path.basename(proteinSequenceFile)}/short_summary.specific.{buscoDataset}.busco_{os.path.basename(proteinSequenceFile)}.txt',
			   'statjson': f'{os.path.dirname(proteinSequenceFile)}/busco_{os.path.basename(proteinSequenceFile)}/short_summary.specific.{buscoDataset}.busco_{os.path.basename(proteinSequenceFile)}.json'}
	protect = [outputs['stattxt'], outputs['statjson']]
	options = {
		'cores': 30,
		'memory': '200g',
		'walltime': '10:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	busco \\
		--cpu {options['cores']} \\
		--force \\
		--in {proteinSequenceFile} \\
		--mode proteins \\
		--out busco_{os.path.basename(proteinSequenceFile)} \\
		--out_path {os.path.dirname(proteinSequenceFile)} \\
		--download_path {buscoDownloadPath} \\
		--lineage {buscoDownloadPath}/lineages/{buscoDataset} \\
		--tar
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec, executor=Conda(environment))

def merqury(genomeAssemblyFile: str, pacbioHifiReads: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genomeAssemblyFile,
		   	  'reads': pacbioHifiReads}
	outputs = {}
	options = {
		'cores': 60,
		'memory': '100g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/merqury ] || mkdir -p {outputDirectory}/merqury
	
	cd {outputDirectory}/merqury

	kmersize="$( \\
	"$(dirname "$(dirname "$(which merqury.sh)")")"/share/merqury/best_k.sh \\
		$(seqtk size {genomeAssemblyFile} | cut -f 2) \\
	| awk \\
		'BEGIN{{FS=OFS=" "}}
		{{if (NR == 3)
			if (($0 - int($0)) >= 0.5)
			{{print int($0) + 1; exit}}
		else
			{{print int($0); exit}}
		}}' \\
	)"

	meryl count \\
		memory={options['memory']} \\
		threads={options['cores']} \\
		k="$kmersize" \\
		output {os.path.basename(os.path.splitext(pacbioHifiReads)[0])}.meryl \\
		{pacbioHifiReads}

	merqury.sh \\
		{os.path.basename(os.path.splitext(pacbioHifiReads)[0])}.meryl \\
		{genomeAssemblyFile} \\
		{os.path.basename(os.path.splitext(pacbioHifiReads)[0])}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def blobtools_blobdir(genomeAssemblyFile: str, speciesName: str, blastnResultFile: str, diamondResultFile: str, coverageAlignmentFile: str, buscoFullTableFile: str, ncbiTaxdumpDirectory: str = "/faststorage/project/EcoGenetics/databases/NCBI_Taxdump"):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genomeAssemblyFile,
			  'blastn': blastnResultFile,
			  'diamond': diamondResultFile,
			  'coverage': coverageAlignmentFile,
			  'busco': buscoFullTableFile}
	outputs = {}
	options = {
		'cores': 32,
		'memory': '80g',
		'walltime': '12:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {os.path.dirname(genomeAssemblyFile)}/qc/blobtools ] || mkdir -p {os.path.dirname(genomeAssemblyFile)}/qc/blobtools

	blobtools create \\
		--threads {options['cores']} \\
		--key assembly.alias="{speciesAbbreviation(speciesName)}" \\
		--key record_type="scaffold" \\
		--key taxon.name="{speciesName}" \\
		--key taxon.genus="{speciesName.split(sep=" ")[0]}" \\
		--key taxon.species="{speciesName}" \\
		--fasta {genomeAssemblyFile} \\
		--hits {blastnResultFile} \\
		--hits {diamondResultFile} \\
		--taxrule bestsumorder \\
		--taxdump {ncbiTaxdumpDirectory} \\
		--cov {coverageAlignmentFile} \\
		--busco {buscoFullTableFile} \\
		{os.path.dirname(genomeAssemblyFile)}/qc/blobtools/blobtools_{os.path.basename(genomeAssemblyFile)}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def blobtools_blastn(genomeAssemblyFile: str, blastDatabase: str = "/faststorage/project/EcoGenetics/databases/NCBI_BLAST_DB/nt/nt"):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genomeAssemblyFile}
	outputs = {'blast': f'{os.path.dirname(genomeAssemblyFile)}/qc/blastn/{os.path.basename(genomeAssemblyFile)}.blast.out'}
	options = {
		'cores': 32,
		'memory': '20g',
		'walltime': '48:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {os.path.dirname(genomeAssemblyFile)}/qc/blastn ] || mkdir -p {os.path.dirname(genomeAssemblyFile)}/qc/blastn
	
	blastn \\
		-num_threads {options['cores']} \\
		-task megablast \\
		-db {blastDatabase} \\
		-query {genomeAssemblyFile} \\
		-outfmt "6 qseqid staxids bitscore std" \\
		-max_target_seqs 10 \\
		-max_hsps 1 \\
		-evalue 1e-25 \\
		-out {os.path.dirname(genomeAssemblyFile)}/qc/blastn/{os.path.basename(genomeAssemblyFile)}.blast.prog.out
	
	mv {os.path.dirname(genomeAssemblyFile)}/qc/blastn/{os.path.basename(genomeAssemblyFile)}.blast.prog.out {outputs['blast']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def blobtools_diamond(genomeAssemblyFile: str, diamondDatabaseFile: str = "/faststorage/project/EcoGenetics/databases/UniProt/reference_proteomes.dmnd"):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genomeAssemblyFile}
	outputs = {'diamond': f'{os.path.dirname(genomeAssemblyFile)}/qc/diamond/{os.path.basename(genomeAssemblyFile)}.diamond.out'}
	options = {
		'cores': 32,
		'memory': '20g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {os.path.dirname(genomeAssemblyFile)}/qc/diamond ] || mkdir -p {os.path.dirname(genomeAssemblyFile)}/qc/diamond
	
	diamond blastx \\
		--threads {options['cores']} \\
		--db {diamondDatabaseFile} \\
		--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
		--sensitive \\
		--max-target-seqs 1 \\
		--evalue 1e-25 \\
		--query {genomeAssemblyFile} \\
		> {os.path.dirname(genomeAssemblyFile)}/qc/diamond/{os.path.basename(genomeAssemblyFile)}.diamond.prog.out
	
	mv {os.path.dirname(genomeAssemblyFile)}/qc/diamond/{os.path.basename(genomeAssemblyFile)}.diamond.prog.out {outputs['diamond']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def blobtools_coverage(genomeAssemblyFile: str, pacbioHifiReads: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genomeAssemblyFile,
		   	  'hifireads': pacbioHifiReads}
	outputs = {'alignment': f'{os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam',
			   'index': f'{os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam.csi'}
	options = {
		'cores': 32,
		'memory': '100g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {os.path.dirname(genomeAssemblyFile)}/qc/coverage ] || mkdir -p {os.path.dirname(genomeAssemblyFile)}/qc/coverage
	
	minimap2 \\
		-x map-hifi \\
		-t {options['cores']} \\
		-a \\
		{genomeAssemblyFile} \\
		{pacbioHifiReads} \\
	| samtools sort \\
		--threads {options['cores'] - 1} \\
		--output-fmt BAM \\
		-o {os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam
	
	samtools index \\
		--threads {options['cores'] - 1} \\
		--csi \\
		--output {os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam.csi \\
		{os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

########################## Purge duplicates ##########################
def purge_dups_1_map_hifi_to_genome(gemoneAssemblyFile: str, hifiSequenceFile: str, outputDirectory: str, speciesName: str, environment: str, directoryAddition: str | None = None, roundNumber: int = 1):
	"""
	Template: Align PacBio HiFi sequence data to genome.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	if directoryAddition:
		directoryAddition = '_' + directoryAddition
	else:
		directoryAddition = ''
	inputs = {'genome': gemoneAssemblyFile,
		   	  'hifi': hifiSequenceFile}
	outputs = {'paf': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.paf.gz',
			   'stat': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/PB.stat',
			   'cov': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/PB.base.cov'}
	options = {
		'cores': 30,
		'memory': '60g',
		'walltime': '04:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments ] || mkdir -p {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments
	
	minimap2 \\
		-x map-hifi \\
		-t {options['cores']} \\
		{gemoneAssemblyFile} \\
		{hifiSequenceFile} \\
	| gzip \\
		-c \\
		- \\
		> {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.prog.paf.gz
	
	pbcstat \\
		-O {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/ \\
		{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.prog.paf.gz
	
	mv {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.prog.paf.gz {outputs['paf']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def purge_dups_2_map_to_self(genomeAssemblyFile: str, outputDirectory: str, speciesName: str, environment: str, directoryAddition: str | None = None, roundNumber: int = 1):
	"""
	Template: Splits assembly into equal pieces and maps it to itself.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	if directoryAddition:
		directoryAddition = '_' + directoryAddition
	else:
		directoryAddition = ''
	inputs = {'genome': genomeAssemblyFile}
	outputs = {'split': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.split.fasta',
			   'paf': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.self.paf.gz'}
	options = {
		'cores': 30,
		'memory': '60g',
		'walltime': '04:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments ] || mkdir -p {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments
	
	split_fa \\
		{genomeAssemblyFile} \\
		> {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.split.prog.fasta
	
	mv {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.split.prog.fasta {outputs['split']}

	minimap2 \\
		-x asm5 \\
		-t {options['cores']} \\
		-D \\
		-P \\
		{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.split.fasta \\
		{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.split.fasta \\
	| gzip \\
		-c \\
		- \\
		> {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.self.prog.paf.gz
	
	mv {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/alignments/{speciesAbbreviation(speciesName)}.self.prog.paf.gz {outputs['paf']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def purge_dups_3_purge_duplicates(pbStatFile: str, pbBaseCovFile: str, selfAlignmentPaf: str, genomeAssemblyFile: str, outputDirectory: str, speciesName: str, environment: str, directoryAddition: str | None = None, roundNumber: int = 1):
	"""
	Template: Purge haplotigs and opverlaps using :script:`purge_dups`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	if directoryAddition:
		directoryAddition = '_' + directoryAddition
	else:
		directoryAddition = ''
	inputs = {'stat': pbStatFile,
		   	  'cov': pbBaseCovFile,
			  'self': selfAlignmentPaf,
			  'genome': genomeAssemblyFile}
	outputs = {'cutoffs': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/cutoffs',
			   'dups': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/dups.bed',
			   'purged': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/{speciesAbbreviation(speciesName)}.purged.fa',
			   'hap': f'{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/{speciesAbbreviation(speciesName)}.hap.fa'}
	options = {
		'cores': 2,
		'memory': '20g',
		'walltime': '04:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02} ] || mkdir -p {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}
	
	calcuts \\
		{pbStatFile} \\
		> {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/cutoffs.prog
	
	mv {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/cutoffs.prog {outputs['cutoffs']}

	purge_dups \\
		-2 \\
		-T {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/cutoffs \\
		-c {pbBaseCovFile} \\
		{selfAlignmentPaf} \\
		> {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/dups.prog.bed

	mv {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/dups.prog.bed {outputs['dups']}

	get_seqs \\
		-p {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/{speciesAbbreviation(speciesName)}.prog \\
		{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/dups.bed \\
		{genomeAssemblyFile}

	seqkit seq \\
		-j {options['cores']} \\
		-w 60 \\
		-o {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/{speciesAbbreviation(speciesName)}.purged.fa \\
		{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/{speciesAbbreviation(speciesName)}.prog.purged.fa

	rm {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/{speciesAbbreviation(speciesName)}.prog.purged.fa
	mv {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/{speciesAbbreviation(speciesName)}.prog.hap.fa {outputs['hap']}
	
	awk \\
		'BEGIN{{
			OFS = "\\t"
		}}
		{{
			if ($0 ~ /^>/)
			{{
				if (sequence_length)
				{{
					print sequence_name, sequence_length
				}}
				total_length += sequence_length
				sequence_name = $0
				sequence_length = 0
				sequence_number += 1
				next
			}}
			sequence_length += length($0)
		}}
		END{{
			if (sequence_length)
			{{
				print sequence_name, sequence_length
			}}
			print sequence_number"_sequences", total_length + sequence_length
		}}' \\
		{outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/{speciesAbbreviation(speciesName)}.purged.fa \\
		> {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/sequences.tsv

	last_line=($(tail -n 1 {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/sequences.tsv))

	mv {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/sequences.tsv {outputDirectory}/purge_dups{directoryAddition}/{roundNumber:02}/sequences_"${{last_line[0]%_*}}"_"${{last_line[1]}}".tsv

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))