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
	outputs = {'stattxt': f'{outputDirectory}/busco_{buscoLineage}/run_{buscoLineage}/short_summary.txt',
			   'statjson': f'{outputDirectory}/busco_{buscoLineage}/run_{buscoLineage}/short_summary.json',
			   'fulltable': f'{outputDirectory}/busco_{buscoLineage}/run_{buscoLineage}/full_table.tsv'}
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
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}

	cd {outputDirectory}

	busco \\
		--cpu {options['cores']} \\
		--force \\
		--metaeuk \\
		--in {genomeAssemblyFile} \\
		--mode genome \\
		--out busco_{buscoLineage} \\
		--out_path {outputDirectory} \\
		--download_path {buscoDownloadPath} \\
		--lineage_dataset {buscoDownloadPath}/lineages/{buscoLineage} \\
		--offline

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_extract_busco_genes(buscoFullTableTsv: str, outputPrefix: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'tsv': buscoFullTableTsv}
	outputs = {'fasta': f'{os.path.dirname(os.path.dirname(buscoFullTableTsv))}/{outputPrefix}.buscoGenes.fasta'}
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
	
	btk pipeline extract-busco-genes \\
		--busco {os.path.dirname(buscoFullTableTsv)}/busco_sequences \\
		--out {os.path.dirname(os.path.dirname(buscoFullTableTsv))}/{outputPrefix}.buscoGenes.prog.fasta

	mv {os.path.dirname(os.path.dirname(buscoFullTableTsv))}/{outputPrefix}.buscoGenes.prog.fasta {outputs['fasta']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def diamond_blastp(queryFileFasta: str, outputDirectory: str, environment: str, diamondDatabaseFile: str = "/faststorage/project/EcoGenetics/databases/UniProt/reference_proteomes.dmnd"):
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
	
	diamond blastp \\
		--threads {options['cores']} \\
		--db {diamondDatabaseFile} \\
		--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
		--evalue 1.0e-25 \\
		--max-target-seq 10 \\
		--max-hsps 1 \\
		--query {queryFileFasta} \\
		--out {outputDirectory}/diamond/blastp/{filename}.diamondBlastp.prog.txt
	
	mv {outputDirectory}/diamond/blastp/{filename}.diamondBlastp.prog.txt {outputs['blast']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def diamond_blastx(queryFileFasta: str, buscoTableFull: str, outputDirectory: str, environment: str, diamondDatabaseFile: str = "/faststorage/project/EcoGenetics/databases/UniProt/reference_proteomes.dmnd"):
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
	outputs = {'blast': f'{outputDirectory}/diamond/blastx/{filename}.diamondBlastx.txt'}
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
	
	[ -d {outputDirectory}/diamond/blastx ] || mkdir -p {outputDirectory}/diamond/blastx
	
	btk pipeline chunk-fasta \\
		--chunk 100000 \\
		--overlap 0 \\
		--max-chunks 10 \\
		--min-length 1000 \\
		--in {queryFileFasta} \\
		--busco {buscoTableFull} \\
		--out {outputDirectory}/diamond/blastx/{filename}.chunks.fasta

	diamond blastx \\
		--threads {options['cores']} \\
		--db {diamondDatabaseFile} \\
		--query {outputDirectory}/diamond/blastx/{filename}.chunks.fasta \\
		--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
		--evalue 1.0e-25 \\
		--max-target-seq 10 \\
		--max-hsps 1 \\
		--out {outputDirectory}/diamond/blastx/{filename}.diamondBlastx.chunks.txt \\
		--log
	
	btk pipeline unchunk-blast \\
		--count 10 \\
		--in {outputDirectory}/diamond/blastx/{filename}.diamondBlastx.chunks.txt \\
		--out {outputDirectory}/diamond/blastx/{filename}.diamondBlastx.prog.txt

	mv {outputDirectory}/diamond/blastx/{filename}.diamondBlastx.prog.txt {outputs['blast']}
		
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

def ncbi_blastn(queryFileFasta: str, outputDirectory: str, environment: str, ncbiBlastDatabase: str = "/faststorage/project/EcoGenetics/databases/NCBI_BLAST_DB/nt/nt"):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(queryFileFasta)[0])[0]) if queryFileFasta.endswith('.gz') else os.path.basename(os.path.splitext(queryFileFasta)[0])
	inputs = {}
	outputs = {'blast': f'{outputDirectory}/blast/blastn/{filename}.blastBlastN.txt'}
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
	
	[ -d {outputDirectory}/blast/blastn ] || mkdir -p {outputDirectory}/blast/blastn
	
	btk pipeline chunk-fasta \\
		--chunk 100000 \\
		--overlap 0 \\
		--max-chunks 10 \\
		--min-length 1000 \\
		--in {queryFileFasta} \\
		--busco None \\
		--out {outputDirectory}/blast/blastn/{filename}.chunks.fasta

	blastn \\
		-num_threads {options['cores']} \\
		-db {ncbiBlastDatabase} \\
		-query {outputDirectory}/blast/blastn/{filename}.chunks.fasta \\
		-task megablast \\
		-outfmt '6 qseqid staxids bitscore std' \\
		-max_target_seqs 10 \\
		-max_hsps 1 \\
		-evalue 1.0e-10 \\
		-lcase_masking \\
		-dust '20 64 1' \\
		-out {outputDirectory}/blast/blastn/{filename}.blastBlastN.chunks.txt \\
        2> >( tee {outputDirectory}/blast/blastn/{filename}.blastBlastN.error.log >&2 ) || true

	if [[ -s {outputDirectory}/blast/blastn/{filename}.blastBlastN.error.log ]]; then
        grep -qF 'BLAST Database error: Taxonomy ID(s) not found.Taxonomy ID(s) not found' {outputDirectory}/blast/blastn/{filename}.error.log
    fi
	
	btk pipeline unchunk-blast \\
		--count 10 \\
		--in {outputDirectory}/blast/blastn/{filename}.blastBlastN.chunks.txt \\
		--out {outputDirectory}/blast/blastn/{filename}.blastBlastN.prog.txt

	mv {outputDirectory}/blast/blastn/{filename}.blastBlastN.prog.txt {outputs['blast']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_collate_stats(genomeAssemblyFile: str, buscoTables: list, windowsBedFile: str, windowFreqFile: str, windowMononucFile: str,
							  regionsCoverageFile: str, outputDirectory: str, environment: str, windowSizes: list = [0.1, 0.01, 1, 100000, 1000000],
							  windowstats_input: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/windowstats_input.py'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(genomeAssemblyFile)[0])[0]) if genomeAssemblyFile.endswith('.gz') else os.path.basename(os.path.splitext(genomeAssemblyFile)[0])
	inputs = {'busco': buscoTables,
		   	  'bed': windowsBedFile,
			  'freq': windowFreqFile,
			  'mononuc': windowMononucFile,
			  'cov': regionsCoverageFile}
	outputs = {'stats': [f'{outputDirectory}/windowStats/{filename}.windowStats.{size}.tsv' for size in windowSizes]}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '12:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/windowStats/tmp ] || mkdir -p {outputDirectory}/windowStats/tmp
	
	btk pipeline count-busco-genes \\
		--in {' '.join(buscoTables)} \\
		--mask {windowsBedFile} \\
		--out {outputDirectory}/windowStats/tmp/{filename}.buscoGenes.tsv
	
	python {windowstats_input} \\
		--freq {windowFreqFile} \\
		--mononuc {windowMononucFile} \\
		--depth {regionsCoverageFile} \\
		--countbusco {outputDirectory}/windowStats/tmp/{filename}.buscoGenes.tsv \\
		--out {outputDirectory}/windowStats/tmp/{filename}.windowStatsInput.tsv
	
	windoiwsSizes=({' '.join(windowSizes)})
	for size in "${{windowSizes[@]}}"; do
		btk pipeline window-stats \\
			--window "$size" \\
			--in {outputDirectory}/windowStats/tmp/{filename}.windowStatsInput.tsv \\
			--out {outputDirectory}/windowStats/{filename}.windowStats."$size".prog.tsv
	done
	
	for size in "${{windowSizes[@]}}"; do
		mv {outputDirectory}/windowStats/{filename}.windowStats."$size".prog.tsv {outputDirectory}/windowStats/{filename}.windowStats."$size".tsv
	done

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))

def blobtoolkit_create_blobdir(windowStatsFiles: list, outputDirectory: str, environment: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {}
	outputs = {}
	options = {
		'cores': 20,
		'memory': '20g',
		'walltime': '12:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	blobtools replcae \\
		--threads {options['cores']} \\
		--bedtsvdir {os.path.dirname(windowStatsFiles[0])} \\
		
		--evalue 1.0e-25 --hit-count 10
	
	mv
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment))


def busco_protein(proteinSequenceFile: str, buscoDataset: str, buscoDownloadPath: str = '/faststorage/project/EcoGenetics/databases/BUSCO'):
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
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

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