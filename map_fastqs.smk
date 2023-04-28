'''
maps fastq files to a reference genome, converts them to bam format, sorts them,
and indexes the sorted bam file.
'''
configfile: 'map_fastqs.yaml'
output_folder=config['output_folder']
rule all:
	input:
		sorted_bam=expand(output_folder+'/bam_main_files/{sample}_sorted.bam', sample=config['samples']),
		sorted_bam_index=expand(output_folder+'/bam_main_files/{sample}_sorted.bam.bai', sample=config['samples'])

rule copy_files:
	'''
	copies snakemake and yaml files to the output folder so users and
	collaborators can see how the data was produced.
	'''
	input:
		snakefile='map_fastqs.smk',
		configfile='map_fastqs.yaml'
	output:
		snakefile=output_folder+'/snakemake_parameters/map_fastqs.smk',
		configfile=output_folder+'/snakemake_parameters/map_fastqs.yaml'
	shell:
		'''
		cp {input.snakefile} {output.snakefile}
		cp {input.configfile} output.configfile}
		'''

if len(config['mate2_suffix'])>0:
	print(len(config['mate2_suffix']))
	rule align_paired_sam:
		'''
		the second fastq file is passed in as a parameter because if you have single
		end sequencing data, then the second fastq file is not a file that 'needs'
		to exist in order for the program to run. 2nd file empty quotes becomes
		extra whitespace.
		'''
		input:
			indexed_genome=config['indexed_genome_path'],
			fastq_file1=config['fastq_folder']+'/{sample}'+config['mate1_suffix'],
			fastq_file2=config['fastq_folder']+'/{sample}'+config['mate2_suffix']
		output:
			sample_sam=output_folder+'/sam_main_files/{sample}.sam'
		conda:
			'envs/bwa.yaml'
		shell:
			'''
			bwa mem -o {output.sample_sam} {input.indexed_genome} {input.fastq_file1} {input.fastq_file2}
			'''

else:
	rule align_nonpaired_sam:
		'''
		the second fastq file is passed in as a parameter because if you have single
		end sequencing data, then the second fastq file is not a file that 'needs'
		to exist in order for the program to run. 2nd file empty quotes becomes
		extra whitespace.
		'''
		input:
			indexed_genome=config['indexed_genome_path'],
			fastq_file1=config['fastq_folder']+'/{sample}'+config['mate1_suffix']
		output:
			sample_sam=output_folder+'/sam_main_files/{sample}.sam'
		conda:
			'envs/bwa.yaml'
		shell:
			'''
			bwa mem -o {output.sample_sam} {input.indexed_genome} {input.fastq_file1}
			'''

rule make_bam:
	input:
		sample_sam=output_folder+'/sam_main_files/{sample}.sam'
	output:
		sample_bam=temp(output_folder+'/bam_temp_files/{sample}.bam')
	conda:
		'envs/samtools.yaml'
	shell:
		'samtools view -b -o {output.sample_bam} {input.sample_sam}'

rule sort_bam:
	input:
		sample_bam=output_folder+'/bam_temp_files/{sample}.bam'
	output:
		sorted_bam=output_folder+'/bam_main_files/{sample}_sorted.bam'
	conda:
		'envs/samtools.yaml'
	shell:
		'samtools sort -o {output.sorted_bam} {input.sample_bam}'

rule index_bam:
	input:
		sorted_bam=output_folder+'/bam_main_files/{sample}_sorted.bam'
	output:
		bam_index=output_folder+'/bam_main_files/{sample}_sorted.bam.bai'
	conda:
		'envs/samtools.yaml'
	shell:
		'samtools index {input.sorted_bam}'
