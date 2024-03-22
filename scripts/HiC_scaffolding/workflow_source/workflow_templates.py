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

#--------------------------------------------------------------------------
########################## Hi-C scaffolding YaHS ##########################
#--------------------------------------------------------------------------

def index_reference(reference_genome_file: str, output_directory: str):
	"""
	Template: Index reference genome using both :script:`bwa index` and :script:`samtools faidx`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference': reference_genome_file}
	outputs = {'bwa': [f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.amb',
					   f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.ann',
					   f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.pac',
					   f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.bwt',
					   f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.sa'],
			   'fai': f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.fai'}
	protect = [outputs['bwa'][0], outputs['bwa'][1], outputs['bwa'][2], outputs['bwa'][3], outputs['bwa'][4], outputs['fai']]
	options = {
		'cores': 1,
		'memory': '10g',
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
	
	[ -d {output_directory}/HiC_scaffolding/YaHS/reference ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/reference
	ln -s {reference_genome_file} {output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}

	bwa index \
		-p {output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)} \
		{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}
	
	samtools faidx \
		-o {output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.prog.fai \
		{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}
	
	mv {output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.prog.fai {outputs['fai']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def hic_alignment_to_draft_assembly(hic_seqeuence_files: list, draft_genome_file: str, reference_indices: list, output_directory: str, species_name: str):
	"""
	Template: Align Hi-C data to draft genome using :script:`bwa mem`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hic_reads': hic_seqeuence_files,
		   	  'genome': draft_genome_file,
			  'indices': reference_indices}
	outputs = {'bam': f'{output_directory}/HiC_scaffolding/YaHS/alignment/{species_abbreviation(species_name)}.HiC_to_draft.bam'}
	options = {
		'cores': 36,
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
	
	[ -d {output_directory}/HiC_scaffolding/YaHS/alignment/tmp ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/alignment/tmp
	
	bwa mem \
        -t {options['cores']} \
        -R '@RG\\tID:{species_abbreviation(species_name)}.HiC\\tSM:HiC_to_draft' \
        -S \
        -P \
        -5 \
        -T 0 \
        {draft_genome_file} \
        {hic_seqeuence_files[0]} \
        {hic_seqeuence_files[1]} \
    | samtools sort \
        -@ {int(options['cores']) - 1} \
		-n \
        -O BAM \
        -T {output_directory}/HiC_scaffolding/YaHS/alignment/tmp \
        -o {output_directory}/HiC_scaffolding/YaHS/alignment/{species_abbreviation(species_name)}.HiC_to_draft.prog.bam
	
	mv {output_directory}/HiC_scaffolding/YaHS/alignment/{species_abbreviation(species_name)}.HiC_to_draft.prog.bam {outputs['bam']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mark_duplicates_picard(alignment_bam_file: str, output_directory: str):
	"""
	Template: Mark duplicates in :format:`BAM` alignment file using PICARD :script:`MarkDuplicates`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'bam': alignment_bam_file}
	outputs = {'markdup': f'{os.path.splitext(alignment_bam_file)[0]}.markdup.bam',
			   'metric': f'{os.path.splitext(alignment_bam_file)[0]}.markdup.txt'}
	protect = outputs['metric']
	options = {
		'cores': 1,
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
	
	[ -d {output_directory}/HiC_scaffolding/YaHS/alignment/tmp ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/alignment/tmp
	
	export _JAVA_OPTIONS="-Xmx{options['memory']}"

	picard MarkDuplicates \
		--INPUT {alignment_bam_file} \
		--OUTPUT {os.path.splitext(alignment_bam_file)[0]}.markdup.prog.bam \
		--METRICS_FILE {os.path.splitext(alignment_bam_file)[0]}.markdup.prog.txt \
		--TMP_DIR {output_directory}/HiC_scaffolding/YaHS/alignment/tmp
	
	mv {os.path.splitext(alignment_bam_file)[0]}.markdup.prog.bam {outputs['markdup']}
	mv {os.path.splitext(alignment_bam_file)[0]}.markdup.prog.txt {outputs['metric']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def hic_scaffolding_yahs(draft_genome_file: str, hic_to_draft_bam_file: str, output_directory: str, species_name: str):
	"""
	Template: Genome scaffolding with Hi-C data using :script:`YaHS`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': draft_genome_file,
		   	  'bam': hic_to_draft_bam_file}
	outputs = {'final_agp': f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}_scaffolds_final.agp',
               'final_fa': f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}_scaffolds_final.fa',
               'bin': f'{output_directory}/HiC/YaHS/{species_abbreviation(species_name)}.bin'}
	protect = [outputs['final_agp'], outputs['final_fa'], outputs['bin']]
	options = {
		'cores': 1,
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
	
	[ -d {output_directory}/HiC_scaffolding/YaHS ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS
	
	yahs \
        -r 500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000 \
        -o {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)} \
        {draft_genome_file} \
        {hic_to_draft_bam_file}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def yahs_conversion_manual_curation(hic_bin_file: str, scaffolds_final_agp_file: str, draft_assembly_fai_index_file: str, output_directory: str, species_name: str):
    """
    Template: Convert Hi-C alignment file to format required by juicer_tools.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'bin': hic_bin_file,
              'scaffolds': scaffolds_final_agp_file,
              'index': draft_assembly_fai_index_file}
    outputs = {'jbat': [f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.txt',
                        f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.liftover.agp',
                        f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.assembly',
                        f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.assembly.agp',
                        f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.log']}
    protect = outputs['jbat']
    options = {
        'cores': 1,
        'memory': '30g',
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
    
    [ -d {output_directory}/HiC_scaffolding/YaHS/tmp ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/tmp
    [ -d {output_directory}/HiC_scaffolding/YaHS ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS

    juicer pre \
        -a \
        -o {output_directory}/HiC_sacffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog \
        {hic_bin_file} \
        {scaffolds_final_agp_file} \
        {draft_assembly_fai_index_file} \
        > {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.log 2>&1
    
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.txt {outputs['jbat'][0]}
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.liftover.agp {outputs['jbat'][1]}
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.assembly {outputs['jbat'][2]}
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.assembly.agp {outputs['jbat'][3]}
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.log {outputs['jbat'][4]}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def chrom_sizes(scaffolds_final_fa_file: str, chrom_size_script: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/get_length.py'):
    """
    Template: Create file containing two columns. In column 1 the name of each sequence, in column 2 the length of each sequence.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'scaffolds': scaffolds_final_fa_file}
    outputs = {'sizes': f'{scaffolds_final_fa_file}.chrom_sizes'}
    protect = outputs['sizes']
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '06:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    python {chrom_size_script} \
        {scaffolds_final_fa_file} \
        {scaffolds_final_fa_file}.prog.chrom_sizes
    
    mv {scaffolds_final_fa_file}.prog.chrom_sizes {outputs['sizes']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def contact_matrix_manual_curation(JBAT_text_file: str, JBAT_log_file: str, output_directory: str, species_name: str, juicer_tools: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/juicer_tools.2.20.00.jar'):
    """
    Template: Generate Hi-C contact matrix using :script:`juicer pre`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'jbat': JBAT_text_file,
              'sizes': JBAT_log_file}
    outputs = {'hic': f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.hic'}
    protect = outputs['hic']
    options = {
        'cores': 32,
        'memory': '100g',
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
    
    [ -d {output_directory}/HiC_scaffolding/YaHS/tmp ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/tmp

    java -jar -Xmx{options['memory']} {juicer_tools} pre \
        -t {output_directory}/HiC_scaffolding/YaHS/tmp \
        -j {options['cores']} \
        {JBAT_text_file} \
        {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.hic \
        <(cat {JBAT_log_file} | grep PRE_C_SIZE | awk '{{print $2" "$3}}')

    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.hic {outputs['hic']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

#----------------------------------------------------------------------------
########################## Hi-C scaffolding juicer ##########################
#----------------------------------------------------------------------------

def setup_for_juicer(hic_sequence_files: list, draft_genome_assembly_file: str, output_directory: str, split_size: int = 8000000):
	"""
	Template: Set up directory structure required by the :script:`juicer` pipeline.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hic_reads': hic_sequence_files,
		   	  'genome': draft_genome_assembly_file}
	outputs = {'reference': f'{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}',
			   'hic_reads': [f'{output_directory}/HiC_scaffolding/juicer/fastq/hic_R1.fastq',
							 f'{output_directory}/HiC_scaffolding/juicer/fastq/hic_R2.fastq'],
			   'bwa': [f'{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}.amb',
					   f'{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}.ann',
					   f'{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}.pac',
					   f'{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}.bwt',
					   f'{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}.sa'],
			   'fai': f'{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}.fai',
			   'sizes': f'{output_directory}/HiC_scaffolding/juicer/chrom.sizes'}
	protect = [outputs['reference'], outputs['hic_reads'][0], outputs['hic_reads'][1], outputs['bwa'][0], outputs['bwa'][1], outputs['bwa'][2], outputs['bwa'][3], outputs['bwa'][4], outputs['fai'], outputs['sizes']]
	options = {
		'cores': 1,
		'memory': '20g',
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
	
	[ -d {output_directory}/HiC_scaffolding ] || mkdir -p {output_directory}/HiC_scaffolding
	[ -d {output_directory}/HiC_scaffolding/juicer/aligned ] && rm -rf {output_directory}/HiC_scaffolding/juicer/aligned
	[ -d {output_directory}/HiC_scaffolding/juicer/debug ] && rm -rf {output_directory}/HiC_scaffolding/juicer/debug
	[ -d {output_directory}/HiC_scaffolding/juicer/HIC_tmp ] && rm -rf {output_directory}/HiC_scaffolding/juicer/HIC_tmp
	[ -d {output_directory}/HiC_scaffolding/juicer/fastq ] || mkdir -p {output_directory}/HiC_scaffolding/juicer/fastq
	[ -d {output_directory}/HiC_scaffolding/juicer/references ] || mkdir -p {output_directory}/HiC_scaffolding/juicer/references
	[ -d {output_directory}/HiC_scaffolding/juicer/restriction_sites ] || mkdir -p {output_directory}/HiC_scaffolding/juicer/restriction_sites
	[ -d {output_directory}/HiC_scaffolding/juicer/splits ] && rm -rf {output_directory}/HiC_scaffolding/juicer/splits 
	[ -d {output_directory}/HiC_scaffolding/juicer/splits ] || mkdir -p {output_directory}/HiC_scaffolding/juicer/splits
	
	[ -e {outputs['reference']} ] || ln -s {draft_genome_assembly_file} {outputs['reference']}
	gunzip -c {hic_sequence_files[0]} > {outputs['hic_reads'][0]}
	gunzip -c {hic_sequence_files[1]} > {outputs['hic_reads'][1]}

	bwa index \
		-p {output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)} \
		{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}
	
	samtools faidx \
		-o {output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}.prog.fai \
		{output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}
	
	mv {output_directory}/HiC_scaffolding/juicer/references/{os.path.basename(draft_genome_assembly_file)}.prog.fai {outputs['fai']}

	split \
		-a 3 \
		-l {split_size} \
		-d --additional-suffix=_R1.fastq \
		{output_directory}/HiC_scaffolding/juicer/fastq/hic_R1.fastq \
		{output_directory}/HiC_scaffolding/juicer/splits/x
	
	split \
		-a 3 \
		-l {split_size} \
		-d --additional-suffix=_R2.fastq \
		{output_directory}/HiC_scaffolding/juicer/fastq/hic_R2.fastq \
		{output_directory}/HiC_scaffolding/juicer/splits/x
	
	awk \
		'BEGIN{{OFS = "\\t"}}
		{{print $1, $2}}' \
		{outputs['fai']} \
		> {output_directory}/HiC_scaffolding/juicer/chrom.prog.sizes

	mv {output_directory}/HiC_scaffolding/juicer/chrom.prog.sizes {outputs['sizes']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def juicer_pipeline(draft_genome_assembly_file: str, chrom_sizes_file: str, output_directory: str, species_name: str, juicer: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/juicer-1.6/scripts/juicer.sh'):
	"""
	Template: Starts the :script:`juicer` pipeline.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': draft_genome_assembly_file,
		   	  'sizes': chrom_sizes_file}
	outputs = {'nodups': f'{output_directory}/HiC_scaffolding/juicer/aligned/merged_nodups.txt'}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '01:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	bash {juicer} \
		-d {output_directory}/HiC_scaffolding/juicer \
		-D {os.path.dirname(os.path.dirname(juicer))} \
		-p {chrom_sizes_file} \
		-A EcoGenetics \
		-s none \
		-z {draft_genome_assembly_file} \
		-q short \
		-Q 12:00:00 \
		-l normal \
		-L 24:00:00 \
		-t 40 \
		> {species_abbreviation(species_name)}_juicer.log
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def assembly_3ddna(draft_assembly_file: str, merged_nodups_file: str, working_directory: str, input_size: int = 15000, edit_rounds: int = 3):
	"""
	Template: 3D de novo assembly
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	working_dir = working_directory
	inputs = {'assembly': draft_assembly_file,
		   	  'nodups': merged_nodups_file}
	outputs = {'final_fasta': f'{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.fasta',
			   'final_asm': f'{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.asm',
			   'final_assembly': f'{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.assembly',
			   'final_cprops': f'{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.cprops',
			   'final_hic': f'{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.hic'}
	options = {
		'cores': 40,
		'memory': '80g',
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
	
	[ -d {working_directory}/tmp ] || mkdir -p {working_directory}/tmp
	
	export _JAVA_OPTIONS="-Djava.io.tmpdir={working_directory}/tmp -Xmx{options['memory']}"

	3d-dna \
		--rounds {edit_rounds} \
		--input {input_size} \
		{draft_assembly_file} \
		{merged_nodups_file}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)

def finalize_3ddna(reviewed_assembly_file: str, draft_assembly_file: str, merged_nodups_file: str, number_of_chromosomes: int, working_directory: str, post_review: str = '/home/jepe/miniconda3/envs/assembly/share/3d-dna/run-asm-pipeline-post-review.sh'):
	"""
	Template: Finalizes draft assembly using reviewed :format:`assembly` file from :script:`JuiceBox`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	working_dir = working_directory
	inputs = {'asm': reviewed_assembly_file,
		   	  'assembly': draft_assembly_file,
			  'nodups': merged_nodups_file}
	outputs = {'fasta': f'{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}_HiC.fasta',
			   'asm': f'{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}_HiC.assembly',
			   'final_fasta': f'{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.nodebris.fasta'}
	options = {
		'cores': 40,
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
	
	bash {post_review} \
		--sort-output \
		-s finalize \
		-r {reviewed_assembly_file} \
		{draft_assembly_file} \
		{merged_nodups_file}

	seqkit range \
		-j {options['cores']} \
		-r 1:{number_of_chromosomes} \
		-o {working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.nodebris.prog.fasta \
		{working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}_HiC.fasta
	
	mv {working_directory}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.nodebris.prog.fasta {outputs['final_fasta']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)