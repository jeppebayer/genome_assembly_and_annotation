#!/bin/env python3
from gwf import AnonymousTarget
from gwf.executors import Conda
import os, yaml

############################## Functions ##############################

def speciesAbbreviation(speciesName: str) -> str:
	"""Creates species abbreviation from species name.

	:param str speciesName:
		Species name written as *genus* *species*"""
	genus, species = speciesName.replace(' ', '_').split('_')
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

############################## Targets ##############################

def genome_size(genomeAssemblyFile: str, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(genomeAssemblyFile)[0])[0]) if genomeAssemblyFile.endswith('.gz') else os.path.basename(os.path.splitext(genomeAssemblyFile)[0])
	inputs = {'genome': genomeAssemblyFile}
	outputs = {'stats': f'{outputDirectory}/genome_assembly/{filename}.summaryStats.tsv'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '00:30:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/genome_assembly ] || mkdir -p {outputDirectory}/genome_assembly
	
	ln \\
		-s \\
		-f \\
		{genomeAssemblyFile} \\
		{outputDirectory}/genome_assembly/{os.path.basename(genomeAssemblyFile)}

	awk \\
		-v byteSize="$( \\
			wc \\
				--bytes \\
				< {genomeAssemblyFile})" \\
		-v filename={genomeAssemblyFile} \\
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
			array["assemblyFile"] = filename
			if ((byteSize / 1000000) % int(byteSize / 1000000) >= 0.5)
			{{
				array["fileSize(MB)"] = int(byteSize / 1000000) + 1
			}}
			else
			{{
				array["fileSize(MB)"] = int(byteSize / 1000000)
			}}
			if (($2 / 1000000) % int($2 / 1000000) >= 0.5)
			{{
				array["genomeSize(Mb)"] = int($2 / 1000000) + 1
			}}
			else
			{{
				array["genomeSize(Mb)"] = int($2 / 1000000)
			}}
			array["nScaffolds"] = $1
		}}
		END{{
			for (i in array)
			{{
				print i, array[i]
			}}
		}}' \\
		<(seqtk size \\
			{genomeAssemblyFile}) \\
		> {outputDirectory}/genome_assembly/{filename}.summaryStats.prog.tsv
	
	mv {outputDirectory}/genome_assembly/{filename}.summaryStats.prog.tsv {outputs['stats']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def taxon_info(querySpecies: str, taxdumpPath: str, requestedBuscos: str, assemblyFile: str, outputDirectory: str, environment: str, taxonInfo: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/taxon_info/taxonInfo.py'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': assemblyFile}
	outputs = {'info': f'{outputDirectory}/{speciesAbbreviation(querySpecies)}.info.yml'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '00:30:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	python {taxonInfo} \\
		'{querySpecies}' \\
		{taxdumpPath} \\
		{requestedBuscos} \\
		{assemblyFile} \\
		{outputDirectory}/{speciesAbbreviation(querySpecies)}.prog
	
	mv {outputDirectory}/{speciesAbbreviation(querySpecies)}.prog.info.yml {outputs['info']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def fasta_windows(genomeAssemblyFile: str, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(genomeAssemblyFile)[0])[0]) if genomeAssemblyFile.endswith('.gz') else os.path.basename(os.path.splitext(genomeAssemblyFile)[0])
	inputs = {'assembly': genomeAssemblyFile}
	outputs = {'freq': f'{outputDirectory}/fw_out/{filename}.freq_windows.tsv',
			   'mono': f'{outputDirectory}/fw_out/{filename}.mononuc_windows.tsv',
			   'di': f'{outputDirectory}/fw_out/{filename}.dinuc_windows.tsv',
			   'tri': f'{outputDirectory}/fw_out/{filename}.trinuc_windows.tsv',
			   'tetra': f'{outputDirectory}/fw_out/{filename}.tetranuc_windows.tsv',
			   'bed': f'{outputDirectory}/fw_out/{filename}.1k.bed'}
	options = {
		'cores': 1,
		'memory': '20g',
		'walltime': '02:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	cd {outputDirectory}

	fasta_windows \\
		--fasta {genomeAssemblyFile} \\
		--output {filename}.prog
	
	awk \\
		'BEGIN {{
			FS = OFS = "\\t"
		}}
		{{
			if (NR == 1)
			{{
				next
			}}
			print $1, $2, $3
		}}' \\
		{outputDirectory}/fw_out/{filename}.prog_mononuc_windows.tsv \\
		> {outputDirectory}/fw_out/{filename}.1k.prog.bed

	mv {outputDirectory}/fw_out/{filename}.prog_freq_windows.tsv {outputs['freq']}
	mv {outputDirectory}/fw_out/{filename}.prog_mononuc_windows.tsv {outputs['mono']}
	mv {outputDirectory}/fw_out/{filename}.prog_dinuc_windows.tsv {outputs['di']}
	mv {outputDirectory}/fw_out/{filename}.prog_trinuc_windows.tsv {outputs['tri']}
	mv {outputDirectory}/fw_out/{filename}.prog_tetranuc_windows.tsv {outputs['tetra']}
	mv {outputDirectory}/fw_out/{filename}.1k.prog.bed {outputs['bed']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_alignment(genomeAssemblyFile: str, pacbioHifiReads: list, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(genomeAssemblyFile)[0])[0]) if genomeAssemblyFile.endswith('.gz') else os.path.basename(os.path.splitext(genomeAssemblyFile)[0])
	inputs = {'assembly': genomeAssemblyFile,
		   	  'hifireads': pacbioHifiReads}
	outputs = {'alignment': f'{outputDirectory}/alignment/{filename}.hifireads.bam',
			   'index': f'{outputDirectory}/alignment/{filename}.hifireads.bam.csi'}
	options = {
		'cores': 32,
		'memory': '100g',
		'walltime': '24:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/alignment/tmp ] || mkdir -p {outputDirectory}/alignment/tmp
	
	minimap2 \\
		-t {options['cores']} \\
		-x map-hifi \\
		-a \\
		--cs=short \\
		{genomeAssemblyFile} \\
		{' '.join(pacbioHifiReads)} \\
	| samtools sort \\
		-@ {options['cores'] - 1} \\
		-O BAM \\
		-T {outputDirectory}/alignment/tmp \\
		-o {outputDirectory}/alignment/{filename}.hifireads.prog.bam \\
		-
	
	samtools index \\
		--threads {options['cores'] - 1} \\
		--csi \\
		{outputDirectory}/alignment/{filename}.hifireads.prog.bam
	
	mv {outputDirectory}/alignment/{filename}.hifireads.prog.bam {outputs['alignment']}
	mv {outputDirectory}/alignment/{filename}.hifireads.prog.bam.csi {outputs['index']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def samtools_flagstat(alignmnetBamFile: str, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(alignmnetBamFile)[0])[0]) if alignmnetBamFile.endswith('.gz') else os.path.basename(os.path.splitext(alignmnetBamFile)[0])
	inputs = {'bam': alignmnetBamFile}
	outputs = {'flagstat': f'{outputDirectory}/alignment/{filename}.flagstat.tsv'}
	options = {
		'cores': 20,
		'memory': '20g',
		'walltime': '02:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/alignment ] || mkdir -p {outputDirectory}/alignment
	
	samtools flagstat \\
		--threads {options['cores'] - 1} \\
		--output-fmt tsv \\
		{alignmnetBamFile} \\
		> {outputDirectory}/alignment/{filename}.flagstat.prog.tsv
	
	mv {outputDirectory}/alignment/{filename}.flagstat.prog.tsv {outputs['flagstat']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_coverage(hifiToAssemblyBam: str, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(hifiToAssemblyBam)[0])[0]) if hifiToAssemblyBam.endswith('.gz') else os.path.basename(os.path.splitext(hifiToAssemblyBam)[0])
	inputs = {'bam': hifiToAssemblyBam}
	outputs = {'bed': f'{outputDirectory}/coverage/{filename}.regionsCoverage.bed'}
	options = {
		'cores': 1,
		'memory': '30g',
		'walltime': '08:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/coverage ] || mkdir -p {outputDirectory}/coverage
	
	blobtk depth \\
		--bin-size 1000 \\
		--bam {hifiToAssemblyBam} \\
		--bed {outputDirectory}/coverage/{filename}.regionsCoverage.prog.bed
	
	mv {outputDirectory}/coverage/{filename}.regionsCoverage.prog.bed {outputs['bed']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def busco_genome(genomeAssemblyFile: str, buscoLineage: str, buscoDownloadPath: str, outputDirectory: str, environment: str):
	"""
	Template: Runs BUSCO analysis on genome assembly.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': genomeAssemblyFile}
	outputs = {'stattxt': f'{outputDirectory}/busco/{buscoLineage}/run_{buscoLineage}/short_summary.txt',
			   'statjson': f'{outputDirectory}/busco/{buscoLineage}/run_{buscoLineage}/short_summary.json',
			   'fulltable': f'{outputDirectory}/busco/{buscoLineage}/run_{buscoLineage}/full_table.tsv'}
	protect = [outputs['stattxt'], outputs['statjson']]
	options = {
		'cores': 30,
		'memory': '200g',
		'walltime': '18:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/busco ] || mkdir -p {outputDirectory}/busco

	cd {outputDirectory}/busco

	busco \\
		--cpu {options['cores']} \\
		--force \\
		--metaeuk \\
		--in {genomeAssemblyFile} \\
		--mode genome \\
		--out {buscoLineage} \\
		--out_path {outputDirectory}/busco \\
		--download_path {buscoDownloadPath} \\
		--lineage_dataset {buscoDownloadPath}/lineages/{buscoLineage} \\
		--offline

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_extract_busco_genes(buscoFullTableTsv: str, outputPrefix: str, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'tsv': buscoFullTableTsv}
	outputs = {'fasta': f'{outputDirectory}/busco_genes/{outputPrefix}_buscoGenes.fasta'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '01:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/busco_genes ] || mkdir -p {outputDirectory}/busco_genes

	btk pipeline extract-busco-genes \\
		--busco {os.path.dirname(buscoFullTableTsv)}/busco_sequences \\
		--out {outputDirectory}/busco_genes/{outputPrefix}_buscoGenes.prog.fasta

	mv {outputDirectory}/busco_genes/{outputPrefix}_buscoGenes.prog.fasta {outputs['fasta']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def diamond_blastp(queryFileFasta: str, excludeTaxon: str, outputDirectory: str, environment: str, diamondDatabaseFile: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(queryFileFasta)[0])[0]) if queryFileFasta.endswith('.gz') else os.path.basename(os.path.splitext(queryFileFasta)[0])
	inputs = {'query': queryFileFasta}
	outputs = {'blast': f'{outputDirectory}/diamond/blastp/{filename}.diamondBlastp.txt'}
	options = {
		'cores': 30,
		'memory': '20g',
		'walltime': '24:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/diamond/blastp ] || mkdir -p {outputDirectory}/diamond/blastp
	
	cd {outputDirectory}/diamond/blastp

	diamond blastp \\
		--threads {options['cores']} \\
		--db {diamondDatabaseFile} \\
		--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
		--evalue 1.0e-25 \\
		--max-target-seqs 10 \\
		--max-hsps 1 \\
		--taxon-exclude {excludeTaxon} \\
		--query {queryFileFasta} \\
		--out {outputDirectory}/diamond/blastp/{filename}.diamondBlastp.prog.txt
	
	mv {outputDirectory}/diamond/blastp/{filename}.diamondBlastp.prog.txt {outputs['blast']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def diamond_blastx(queryFileFasta: str, buscoTableFull: str, buscoLineage: str, excludeTaxon: str, outputDirectory: str, environment: str, diamondDatabaseFile: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(queryFileFasta)[0])[0]) if queryFileFasta.endswith('.gz') else os.path.basename(os.path.splitext(queryFileFasta)[0])
	inputs = {'query': queryFileFasta,
		   	  'busco': buscoTableFull}
	outputs = {'blast': f'{outputDirectory}/diamond/blastx/{filename}.{buscoLineage}.diamondBlastx.txt'}
	options = {
		'cores': 30,
		'memory': '60g',
		'walltime': '24:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/diamond/blastx ] || mkdir -p {outputDirectory}/diamond/blastx
	
	cd {outputDirectory}/diamond/blastx
	
	btk pipeline chunk-fasta \\
		--chunk 100000 \\
		--overlap 0 \\
		--max-chunks 10 \\
		--min-length 1000 \\
		--in {queryFileFasta} \\
		--busco {buscoTableFull} \\
		--out {outputDirectory}/diamond/blastx/{filename}.{buscoLineage}.chunks.fasta

	diamond blastx \\
		--threads {options['cores']} \\
		--db {diamondDatabaseFile} \\
		--query {outputDirectory}/diamond/blastx/{filename}.{buscoLineage}.chunks.fasta \\
		--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
		--evalue 1.0e-25 \\
		--max-target-seqs 10 \\
		--max-hsps 1 \\
		--taxon-exclude {excludeTaxon} \\
		--out {outputDirectory}/diamond/blastx/{filename}.{buscoLineage}.diamondBlastx.chunks.txt \\
		--log
	
	btk pipeline unchunk-blast \\
		--count 10 \\
		--in {outputDirectory}/diamond/blastx/{filename}.{buscoLineage}.diamondBlastx.chunks.txt \\
		--out {outputDirectory}/diamond/blastx/{filename}.{buscoLineage}.diamondBlastx.prog.txt

	mv {outputDirectory}/diamond/blastx/{filename}.{buscoLineage}.diamondBlastx.prog.txt {outputs['blast']}
		
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def diamond_blastx_no_hits(queryFileFasta: str, blastxResults: str, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(queryFileFasta)[0])[0]) if queryFileFasta.endswith('.gz') else os.path.basename(os.path.splitext(queryFileFasta)[0])
	inputs = {'query': queryFileFasta,
		   	  'blastx': blastxResults}
	outputs = {'txt': f'{outputDirectory}/diamond/blastx/{filename}.nohit.txt',
			   'fasta': f'{outputDirectory}/diamond/blastx/{filename}.nohit.fa'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '06:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/diamond/blastx ] || mkdir -p {outputDirectory}/diamond/blastx
	
	grep \\
		'>' \\
		{queryFileFasta} \\
	| grep \\
		-v \\
		-w \\
		-f <(\\
			awk \\
				-v evalue="1.0e-25" \\
				'{{
					if($14 < evalue)
					{{
						print $1
					}}
				}}' \\
				{blastxResults} \\
			| sort \\
			| uniq \\
			) \\
	| awk \\
		'{{
			print substr($1, 2)
		}}' \\
		> {outputDirectory}/diamond/blastx/{filename}.nohit.txt

	seqtk subseq \\
		{queryFileFasta} \\
		{outputDirectory}/diamond/blastx/{filename}.nohit.txt \\
		> {outputDirectory}/diamond/blastx/{filename}.nohit.prog.fa
	
	mv {outputDirectory}/diamond/blastx/{filename}.nohit.prog.fa {outputs['fasta']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def ncbi_blastn(queryFileFasta: str, outputDirectory: str, excludeTaxon: str, environment: str, ncbiBlastDatabase: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(queryFileFasta)[0])[0]) if queryFileFasta.endswith('.gz') else os.path.basename(os.path.splitext(queryFileFasta)[0])
	tag = ''
	if excludeTaxon:
		tag = f'.exclude{excludeTaxon}'
		excludeTaxon = f'-negative_taxids {excludeTaxon}'
	inputs = {'query': queryFileFasta}
	outputs = {'blast': f'{outputDirectory}/blast/blastn/{filename}{tag}.blastBlastN.txt'}
	options = {
		'cores': 40,
		'memory': '20g',
		'walltime': '72:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/blast/blastn ] || mkdir -p {outputDirectory}/blast/blastn
	
	export BLASTDB="{ncbiBlastDatabase}"

	btk pipeline chunk-fasta \\
		--chunk 100000 \\
		--overlap 0 \\
		--max-chunks 10 \\
		--min-length 1000 \\
		--in {queryFileFasta} \\
		--busco None \\
		--out {outputDirectory}/blast/blastn/{filename}{tag}.chunks.fasta

	blastn \\
		-num_threads {options['cores']} \\
		-db {ncbiBlastDatabase}/nt \\
		-query {outputDirectory}/blast/blastn/{filename}{tag}.chunks.fasta \\
		-task megablast \\
		-outfmt '6 qseqid staxids bitscore std' \\
		-max_target_seqs 10 \\
		-max_hsps 1 \\
		-evalue 1.0e-10 \\
		-lcase_masking \\
		-dust '20 64 1' \\
		{excludeTaxon} \\
		-out {outputDirectory}/blast/blastn/{filename}{tag}.blastBlastN.chunks.txt \\
		2> >( tee {outputDirectory}/blast/blastn/{filename}{tag}.blastBlastN.error.log >&2 ) || true

	if [[ -s {outputDirectory}/blast/blastn/{filename}{tag}.blastBlastN.error.log ]]; then
		grep -qF 'BLAST Database error: Taxonomy ID(s) not found.Taxonomy ID(s) not found' {outputDirectory}/blast/blastn/{filename}{tag}.error.log
	fi
	
	btk pipeline unchunk-blast \\
		--count 10 \\
		--in {outputDirectory}/blast/blastn/{filename}{tag}.blastBlastN.chunks.txt \\
		--out {outputDirectory}/blast/blastn/{filename}{tag}.blastBlastN.prog.txt

	mv {outputDirectory}/blast/blastn/{filename}{tag}.blastBlastN.prog.txt {outputs['blast']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_count_busco_genes(buscoTablesFull: list, windowsBedFile: str, filename: str, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'tables': buscoTablesFull,
		   	  'bed': windowsBedFile}
	outputs = {'count': f'{outputDirectory}/busco_count/{filename}.buscoGenes.tsv'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '00:30:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/busco_count ] || mkdir -p {outputDirectory}/busco_count
	
	btk pipeline count-busco-genes \\
		--in {' --in '.join(buscoTablesFull)} \\
		--mask {windowsBedFile} \\
		--out {outputDirectory}/busco_count/{filename}.prog.buscoGenes.tsv
	
	mv {outputDirectory}/busco_count/{filename}.prog.buscoGenes.tsv {outputs['count']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_windowstats_input(countBuscoGenesFile: str, windowFreqFile: str, windowMononucFile: str, regionsCoverageFile: str,
							  	  filename: str, outputDirectory: str, environment: str, windowstats_input: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/windowstatsInput.py'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genes': countBuscoGenesFile,
			  'freq': windowFreqFile,
			  'mononuc': windowMononucFile,
			  'cov': regionsCoverageFile}
	outputs = {'input': f'{outputDirectory}/windowStats/tmp/{filename}.windowStatsInput.tsv'}
	options = {
		'cores': 1,
		'memory': '20g',
		'walltime': '00:30:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/windowStats/tmp ] || mkdir -p {outputDirectory}/windowStats/tmp
	
	python {windowstats_input} \\
		--freq {windowFreqFile} \\
		--mononuc {windowMononucFile} \\
		--depth {regionsCoverageFile} \\
		--countbuscos {countBuscoGenesFile} \\
		--prefix {outputDirectory}/windowStats/tmp/{filename}.windowStatsInput.prog
	
	mv {outputDirectory}/windowStats/tmp/{filename}.windowStatsInput.prog.tsv {outputs['input']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_windowstats(windowstatsInputFile: str, filename: str, outputDirectory: str, environment: str, windowSizes: list = ['0.1', '0.01', '1', '100000', '1000000']):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'input': windowstatsInputFile}
	outputs = {'stats': f'{outputDirectory}/windowStats/{filename}.windowStats.tsv'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '02:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/windowStats ] || mkdir -p {outputDirectory}/windowStats
	
	btk pipeline window-stats \\
		--window {' --window '.join(windowSizes)} \\
		--in {outputDirectory}/windowStats/tmp/{filename}.windowStatsInput.tsv \\
		--out {outputDirectory}/windowStats/{filename}.windowStats.tsv
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_create_blobdir(genomeAssemblyFile: str, windowstatsFile: str, infoFile: str, ncbiTaxdump: str, diamondBlastpResults: list, buscoTablesFull: list, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(genomeAssemblyFile)[0])[0]) if genomeAssemblyFile.endswith('.gz') else os.path.basename(os.path.splitext(genomeAssemblyFile)[0])
	inputs = {'windowstats': windowstatsFile,
		   	  'info': infoFile,
			  'blastp': diamondBlastpResults,
			  'busco': buscoTablesFull}
	outputs = {'meta': f'{outputDirectory}/blobDir_{filename}/meta.json'}
	options = {
		'cores': 20,
		'memory': '20g',
		'walltime': '02:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/blobDir_{filename} ] || mkdir -p {outputDirectory}/blobDir_{filename}

	blobtools replace \\
		--threads {options['cores']} \\
		--bedtsvdir {os.path.dirname(windowstatsFile)} \\
		--meta {infoFile} \\
		--taxdump {ncbiTaxdump} \\
		--taxrule buscogenes \\
		--hits {' --hits '.join(diamondBlastpResults)} \\
		--busco {' --busco '.join(buscoTablesFull)} \\
		--evalue 1.0e-25 \\
		--hit-count 10 \\
		{outputDirectory}/blobDir_{filename}

	touch {outputDirectory}/blobDir_{filename}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_update_blobdir(diamondBlastxResults: list, ncbiBlastnResults: list, ncbiTaxdump: str, blobdirMeta: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'blastx': diamondBlastxResults,
		   	  'blastn': ncbiBlastnResults,
			  'meta': blobdirMeta}
	outputs = {'buscoregions': f'{os.path.dirname(blobdirMeta)}/buscoregions_kingdom.json'}
	options = {
		'cores': 20,
		'memory': '20g',
		'walltime': '02:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	blobtools replace \\
		--threads {options['cores']} \\
		--taxdump {ncbiTaxdump} \\
		--taxrule bestdistorder=buscoregions \\
		--hits {' --hits '.join(diamondBlastxResults)} \\
		--hits {' --hits '.join(ncbiBlastnResults)} \\
		--evalue 1.0e-25 \\
		--hit-count 10 \\
		{os.path.dirname(blobdirMeta)}
	
	touch {os.path.dirname(blobdirMeta)}/buscoregions_*
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_summary(blobdirFile: str, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'blobdir': blobdirFile}
	outputs = {'summary': f'{os.path.dirname(blobdirFile)}/summary.json'}
	options = {
		'cores': 1,
		'memory': '20g',
		'walltime': '02:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	blobtools filter \\
		--summary {os.path.dirname(blobdirFile)}/summary.prog.json \\
		{os.path.dirname(blobdirFile)}
	
	mv {os.path.dirname(blobdirFile)}/summary.prog.json {outputs['summary']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_images(genomeAssemblyFile: str, blobdirFile: str, outputDirectory: str, environment: str, plots: list = ['blob', 'cumulative', 'snail']):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(genomeAssemblyFile)[0])[0]) if genomeAssemblyFile.endswith('.gz') else os.path.basename(os.path.splitext(genomeAssemblyFile)[0])
	inputs = {'assembly': genomeAssemblyFile,
		   	  'blobdir':blobdirFile}
	outputs = {'plots': [f'{outputDirectory}/plots/{filename}.{plotType}.svg' for plotType in plots]}
	options = {
		'cores': 1,
		'memory': '20g',
		'walltime': '02:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/plots ] || mkdir -p {outputDirectory}/plots
	
	plots=({' '.join(plots)})
	for plotType in "${{plots[@]}}"; do
		if [ "$plotType" == "snail" ]; then
			legendType="default"
		else
			legendType="full"
		fi

		blobtk plot \\
			--view "$plotType" \\
			--blobdir {os.path.dirname(blobdirFile)} \\
			--output {outputDirectory}/plots/{filename}."$plotType".prog.svg \\
			--legend "$legendType"
	done

	for plotType in "${{plots[@]}}"; do
		mv {outputDirectory}/plots/{filename}."$plotType".prog.svg {outputDirectory}/plots/{filename}."$plotType".svg
	done
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def get_versions(versionFiles: list, environment: str, group: str | None = None):
	"""
	Template: Collect version information from all jobs
	
	Template I/O::
	
		inputs = {'versions': list}
		outputs = {'versions': str}
	
	:param list versionFiles:
		List of version files from all jobs.
	"""
	inputs = {'versions': versionFiles}
	outputs = {'versions': f'{os.getcwd()}/versions.yml'}
	options = {
		'cores': 1,
		'memory': '2g',
		'walltime': '00:10:00'
	}
	spec = f"""
	
	cat <<-VERSIONS > {outputs['versions']}
	workflow:
		gwf: $(gwf --version | sed 's/gwf, version //')
	VERSIONS

	cat \\
		{' '.join(versionFiles)} \\
		>> {os.getcwd()}/versions.yml

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

# def busco_protein(proteinSequenceFile: str, buscoDataset: str, buscoDownloadPath: str = '/faststorage/project/EcoGenetics/databases/BUSCO'):
# 	"""
# 	Template: Runs BUSCO analysis on protein sequences from an annotated gene set.
	
# 	Template I/O::
	
# 		inputs = {}
# 		outputs = {}
	
# 	:param
# 	"""
# 	inputs = {'genome': proteinSequenceFile}
# 	outputs = {'stattxt': f'{os.path.dirname(proteinSequenceFile)}/busco_{os.path.basename(proteinSequenceFile)}/short_summary.specific.{buscoDataset}.busco_{os.path.basename(proteinSequenceFile)}.txt',
# 			   'statjson': f'{os.path.dirname(proteinSequenceFile)}/busco_{os.path.basename(proteinSequenceFile)}/short_summary.specific.{buscoDataset}.busco_{os.path.basename(proteinSequenceFile)}.json'}
# 	protect = [outputs['stattxt'], outputs['statjson']]
# 	options = {
# 		'cores': 30,
# 		'memory': '50g',
# 		'walltime': '10:00:00'
# 	}
# 	spec = f"""
# 	# Sources environment
# 	if [ "$USER" == "jepe" ]; then
# 		source /home/"$USER"/.bashrc
# 		source activate assembly
# 	fi
	
# 	echo "START: $(date)"
# 	echo "JobID: $SLURM_JOBID"
	
# 	busco \\
# 		--cpu {options['cores']} \\
# 		--force \\
# 		--in {proteinSequenceFile} \\
# 		--mode proteins \\
# 		--out busco_{os.path.basename(proteinSequenceFile)} \\
# 		--out_path {os.path.dirname(proteinSequenceFile)} \\
# 		--download_path {buscoDownloadPath} \\
# 		--lineage {buscoDownloadPath}/lineages/{buscoDataset} \\
# 		--tar
	
# 	echo "END: $(date)"
# 	echo "$(jobinfo "$SLURM_JOBID")"
# 	"""
# 	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# def merqury(genomeAssemblyFile: str, pacbioHifiReads: str, outputDirectory: str):
# 	"""
# 	Template: template_description
	
# 	Template I/O::
	
# 		inputs = {}
# 		outputs = {}
	
# 	:param
# 	"""
# 	inputs = {'assembly': genomeAssemblyFile,
# 		   	  'reads': pacbioHifiReads}
# 	outputs = {}
# 	options = {
# 		'cores': 60,
# 		'memory': '100g',
# 		'walltime': '24:00:00'
# 	}
# 	spec = f"""
# 	# Sources environment
# 	if [ "$USER" == "jepe" ]; then
# 		source /home/"$USER"/.bashrc
# 		source activate assembly
# 	fi
	
# 	echo "START: $(date)"
# 	echo "JobID: $SLURM_JOBID"
	
# 	[ -d {outputDirectory}/merqury ] || mkdir -p {outputDirectory}/merqury
	
# 	cd {outputDirectory}/merqury

# 	kmersize="$( \\
# 	"$(dirname "$(dirname "$(which merqury.sh)")")"/share/merqury/best_k.sh \\
# 		$(seqtk size {genomeAssemblyFile} | cut -f 2) \\
# 	| awk \\
# 		'BEGIN{{FS=OFS=" "}}
# 		{{if (NR == 3)
# 			if (($0 - int($0)) >= 0.5)
# 			{{print int($0) + 1; exit}}
# 		else
# 			{{print int($0); exit}}
# 		}}' \\
# 	)"

# 	meryl count \\
# 		memory={options['memory']} \\
# 		threads={options['cores']} \\
# 		k="$kmersize" \\
# 		output {os.path.basename(os.path.splitext(pacbioHifiReads)[0])}.meryl \\
# 		{pacbioHifiReads}

# 	merqury.sh \\
# 		{os.path.basename(os.path.splitext(pacbioHifiReads)[0])}.meryl \\
# 		{genomeAssemblyFile} \\
# 		{os.path.basename(os.path.splitext(pacbioHifiReads)[0])}
	
# 	echo "END: $(date)"
# 	echo "$(jobinfo "$SLURM_JOBID")"
# 	"""
# 	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def blobtools_blobdir(genomeAssemblyFile: str, speciesName: str, blastnResultFile: str, diamondResultFile: str, coverageAlignmentFile: str, buscoFullTableFile: str, ncbiTaxdumpDirectory: str = "/faststorage/project/EcoGenetics/databases/NCBI_Taxdump"):
# 	"""
# 	Template: template_description
	
# 	Template I/O::
	
# 		inputs = {}
# 		outputs = {}
	
# 	:param
# 	"""
# 	inputs = {'assembly': genomeAssemblyFile,
# 			  'blastn': blastnResultFile,
# 			  'diamond': diamondResultFile,
# 			  'coverage': coverageAlignmentFile,
# 			  'busco': buscoFullTableFile}
# 	outputs = {}
# 	options = {
# 		'cores': 32,
# 		'memory': '80g',
# 		'walltime': '12:00:00'
# 	}
# 	spec = f"""
# 	# Sources environment
# 	if [ "$USER" == "jepe" ]; then
# 		source /home/"$USER"/.bashrc
# 		source activate assembly
# 	fi
	
# 	echo "START: $(date)"
# 	echo "JobID: $SLURM_JOBID"
	
# 	[ -d {os.path.dirname(genomeAssemblyFile)}/qc/blobtools ] || mkdir -p {os.path.dirname(genomeAssemblyFile)}/qc/blobtools

# 	blobtools create \\
# 		--threads {options['cores']} \\
# 		--key assembly.alias="{speciesAbbreviation(speciesName)}" \\
# 		--key record_type="scaffold" \\
# 		--key taxon.name="{speciesName}" \\
# 		--key taxon.genus="{speciesName.split(sep=" ")[0]}" \\
# 		--key taxon.species="{speciesName}" \\
# 		--fasta {genomeAssemblyFile} \\
# 		--hits {blastnResultFile} \\
# 		--hits {diamondResultFile} \\
# 		--taxrule bestsumorder \\
# 		--taxdump {ncbiTaxdumpDirectory} \\
# 		--cov {coverageAlignmentFile} \\
# 		--busco {buscoFullTableFile} \\
# 		{os.path.dirname(genomeAssemblyFile)}/qc/blobtools/blobtools_{os.path.basename(genomeAssemblyFile)}
	
# 	echo "END: $(date)"
# 	echo "$(jobinfo "$SLURM_JOBID")"
# 	"""
# 	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def blobtools_blastn(genomeAssemblyFile: str, blastDatabase: str = "/faststorage/project/EcoGenetics/databases/NCBI_BLAST_DB/nt/nt"):
# 	"""
# 	Template: template_description
	
# 	Template I/O::
	
# 		inputs = {}
# 		outputs = {}
	
# 	:param
# 	"""
# 	inputs = {'assembly': genomeAssemblyFile}
# 	outputs = {'blast': f'{os.path.dirname(genomeAssemblyFile)}/qc/blastn/{os.path.basename(genomeAssemblyFile)}.blast.out'}
# 	options = {
# 		'cores': 32,
# 		'memory': '20g',
# 		'walltime': '48:00:00'
# 	}
# 	spec = f"""
# 	# Sources environment
# 	if [ "$USER" == "jepe" ]; then
# 		source /home/"$USER"/.bashrc
# 		source activate assembly
# 	fi
	
# 	echo "START: $(date)"
# 	echo "JobID: $SLURM_JOBID"
	
# 	[ -d {os.path.dirname(genomeAssemblyFile)}/qc/blastn ] || mkdir -p {os.path.dirname(genomeAssemblyFile)}/qc/blastn
	
# 	blastn \\
# 		-num_threads {options['cores']} \\
# 		-task megablast \\
# 		-db {blastDatabase} \\
# 		-query {genomeAssemblyFile} \\
# 		-outfmt "6 qseqid staxids bitscore std" \\
# 		-max_target_seqs 10 \\
# 		-max_hsps 1 \\
# 		-evalue 1e-25 \\
# 		-out {os.path.dirname(genomeAssemblyFile)}/qc/blastn/{os.path.basename(genomeAssemblyFile)}.blast.prog.out
	
# 	mv {os.path.dirname(genomeAssemblyFile)}/qc/blastn/{os.path.basename(genomeAssemblyFile)}.blast.prog.out {outputs['blast']}
	
# 	echo "END: $(date)"
# 	echo "$(jobinfo "$SLURM_JOBID")"
# 	"""
# 	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def blobtools_diamond(genomeAssemblyFile: str, diamondDatabaseFile: str = "/faststorage/project/EcoGenetics/databases/UniProt/reference_proteomes.dmnd"):
# 	"""
# 	Template: template_description
	
# 	Template I/O::
	
# 		inputs = {}
# 		outputs = {}
	
# 	:param
# 	"""
# 	inputs = {'assembly': genomeAssemblyFile}
# 	outputs = {'diamond': f'{os.path.dirname(genomeAssemblyFile)}/qc/diamond/{os.path.basename(genomeAssemblyFile)}.diamond.out'}
# 	options = {
# 		'cores': 32,
# 		'memory': '20g',
# 		'walltime': '24:00:00'
# 	}
# 	spec = f"""
# 	# Sources environment
# 	if [ "$USER" == "jepe" ]; then
# 		source /home/"$USER"/.bashrc
# 		source activate assembly
# 	fi
	
# 	echo "START: $(date)"
# 	echo "JobID: $SLURM_JOBID"
	
# 	[ -d {os.path.dirname(genomeAssemblyFile)}/qc/diamond ] || mkdir -p {os.path.dirname(genomeAssemblyFile)}/qc/diamond
	
# 	diamond blastx \\
# 		--threads {options['cores']} \\
# 		--db {diamondDatabaseFile} \\
# 		--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
# 		--sensitive \\
# 		--max-target-seqs 1 \\
# 		--evalue 1e-25 \\
# 		--query {genomeAssemblyFile} \\
# 		> {os.path.dirname(genomeAssemblyFile)}/qc/diamond/{os.path.basename(genomeAssemblyFile)}.diamond.prog.out
	
# 	mv {os.path.dirname(genomeAssemblyFile)}/qc/diamond/{os.path.basename(genomeAssemblyFile)}.diamond.prog.out {outputs['diamond']}
	
# 	echo "END: $(date)"
# 	echo "$(jobinfo "$SLURM_JOBID")"
# 	"""
# 	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def blobtools_coverage(genomeAssemblyFile: str, pacbioHifiReads: str):
# 	"""
# 	Template: template_description
	
# 	Template I/O::
	
# 		inputs = {}
# 		outputs = {}
	
# 	:param
# 	"""
# 	inputs = {'assembly': genomeAssemblyFile,
# 		   	  'hifireads': pacbioHifiReads}
# 	outputs = {'alignment': f'{os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam',
# 			   'index': f'{os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam.csi'}
# 	options = {
# 		'cores': 32,
# 		'memory': '100g',
# 		'walltime': '24:00:00'
# 	}
# 	spec = f"""
# 	# Sources environment
# 	if [ "$USER" == "jepe" ]; then
# 		source /home/"$USER"/.bashrc
# 		source activate assembly
# 	fi
	
# 	echo "START: $(date)"
# 	echo "JobID: $SLURM_JOBID"
	
# 	[ -d {os.path.dirname(genomeAssemblyFile)}/qc/coverage ] || mkdir -p {os.path.dirname(genomeAssemblyFile)}/qc/coverage
	
# 	minimap2 \\
# 		-x map-hifi \\
# 		-t {options['cores']} \\
# 		-a \\
# 		{genomeAssemblyFile} \\
# 		{pacbioHifiReads} \\
# 	| samtools sort \\
# 		--threads {options['cores'] - 1} \\
# 		--output-fmt BAM \\
# 		-o {os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam
	
# 	samtools index \\
# 		--threads {options['cores'] - 1} \\
# 		--csi \\
# 		--output {os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam.csi \\
# 		{os.path.dirname(genomeAssemblyFile)}/qc/coverage/{os.path.basename(genomeAssemblyFile)}.hifireads.bam
	
# 	echo "END: $(date)"
# 	echo "$(jobinfo "$SLURM_JOBID")"
# 	"""
# 	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)