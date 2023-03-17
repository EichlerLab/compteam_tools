import pandas as pd

manifest_df = pd.read_csv('fastq_manifest.tab', sep='\t', index_col='sample', header=0, dtype=str)


def findCram(wildcards):
	return manifest_df.at[wildcards.sample, 'cram']


if 'REF' in config:
	REF_ARGS=f" -T {config['REF']} "
else:
	REF_ARGS=""

localrules: all, fofn

wildcard_constraints:
	sample='|'.join(manifest_df.index)

rule all:
	input:
		expand('fastq/fofn/{sample}_KG.fastq.fofn', sample=manifest_df.index)

rule convert:
	input:
		cram = findCram
	output:
		fastq_one = temp('fastq/ILLUMINA_{sample}_1.fastq'),
		fastq_two = temp('fastq/ILLUMINA_{sample}_2.fastq')
	threads: 8
	resources:
		mem = 12,
		hrs = 48
	shell:
		'''
		module load samtools/1.14
		samtools view -b {REF_ARGS} -@ {threads} {input.cram} | samtools sort -n -T {resources.tmpdir} -@ {threads} | samtools fastq -@ {threads} -1 {output.fastq_one} -2 {output.fastq_two} -0 /dev/null -s /dev/null -n
		'''

rule bgzip:
	input:
		fastq = 'fastq/ILLUMINA_{sample}_{pair}.fastq'
	output:
		gzipped = 'fastq/ILLUMINA_{sample}_{pair}.fastq.gz'
	threads: 1
	resources:
		mem = 12,
		hrs = 24
	shell:
		'''
		module load htslib/1.12
		bgzip -c {input.fastq} > {output.gzipped}
		'''

rule fofn:
	input:
		fastq = expand('fastq/ILLUMINA_{{sample}}_{pair}.fastq.gz', pair=['1','2'])
	output:
		fofn = 'fastq/fofn/{sample}_KG.fastq.fofn'
	shell:
		'''
		readlink -f {input.fastq} > {output.fofn}
		'''
