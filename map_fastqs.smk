'''
maps fastq files to a reference genome, converts them to bam format, sorts them,
and indexes the sorted bam file.
'''
configfile: 'map_fastqs.yaml'
output_folder=config['output_folder']
rule all:
	input:
		sorted_bam=output_folder+'/{sample}_sorted.bam', sample=config['sample_name'],
		sorted_bam_index=output_folder+'/{sample}_sorted.bam.bai', sample=config['sample_name']

rule align_sam:
	input:
		indexed_genome=config['bt2_genome_path'],
		fastq_file1=config['fastq_path_mate1'],
		fastq_file2=config['fastq_path_mate2']
	output:
		sample_sam=output_folder+'/sam_main_files/{sample}.sam'
	conda:
		'envs/bwa.yaml'
	shell:
		'bwa mem -o {output.sample_sam} {input.indexed_genome} {input.fastq_file1} {input.fastq_file2}'

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
