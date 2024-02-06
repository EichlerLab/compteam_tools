import pandas as pd
from pybedtools import BedTool
import numpy as np
import networkx as nx
import more_itertools as mit


configfile: "config/saffire_gen.yaml"


PARTS = config.get("PARTS", 15)
# MINIMAP_PARAMS = config.get('MINIMAP_PARAMS', '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5')
MINIMAP_PARAMS = config.get(
    "MINIMAP_PARAMS",
    "-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5",
)
MANIFEST = config.get("MANIFEST", "config/manifest.tab")
SV_SIZE = config.get("SV_SIZE", "30000")
REF = config.get("REF")

manifest_df = pd.read_csv(MANIFEST, sep="\t", index_col=["SAMPLE"])


def find_fasta(wildcards):
    return manifest_df.loc[wildcards.sample]["ASM"]


def find_contigs(wildcards):
    return gather.split(
        "tmp/{sample}.{scatteritem}-broken.paf", sample=wildcards.sample
    )


scattergather:
    split=PARTS,


wildcard_constraints:
    sample="|".join(manifest_df.index),


localrules:
    all,


rule all:
    input:
        expand("results/saffire/{sample}/{sample}.saf", sample=manifest_df.index),


rule gather_bam:
    input:
        expand("results/{sample}.bam", sample=manifest_df.index),


rule gather_bed:
    input:
        expand("results/{sample}.bed", sample=manifest_df.index),


rule gather_paf:
    input:
        expand("results/{sample}.paf", sample=manifest_df.index),


rule make_paf:
    input:
        fa=find_fasta,
        ref=REF,
    output:
        paf="results/{sample}.paf",
    threads: 8
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "minimap2/2.24",
    resources:
        mem=12,
        hrs=72,
    params:
        minimap=MINIMAP_PARAMS,
    shell:
        """
        minimap2 -c -t {threads} -K {resources.mem}G --eqx --cs {params.minimap} --secondary=no --eqx -Y {input.ref} {input.fa} > {output.paf}
        """


rule make_sam:
    input:
        fa=find_fasta,
        ref=REF,
    output:
        paf="results/{sample}.bam",
    threads: 8
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "minimap2/2.24",
        "samtools/1.12",
    resources:
        mem=60,
        hrs=120,
    params:
        minimap=MINIMAP_PARAMS,
    shell:
        """
        minimap2 -c -t {threads} -K {resources.mem}G --cs -a {params.minimap} --MD --secondary=no --eqx -Y {input.ref} {input.fa} | samtools view -S -b /dev/stdin | samtools sort -@ {threads} /dev/stdin > {output.paf}
        """


rule make_bed:
    input:
        bam=rules.make_sam.output.paf,
    output:
        bed="results/{sample}.bed",
    threads: 1
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "minimap2/2.24",
        "bedtools/2.29.0",
    resources:
        mem=12,
        hrs=72,
    shell:
        """
        bedtools bamtobed -i {input.bam} | bedtools sort -i - | cut -f 1,2,3,4,5 > {output.bed}
        """


rule split_paf:
    input:
        paf=rules.make_paf.output.paf,
    output:
        flag=temp(scatter.split("tmp/{{sample}}.{scatteritem}.paf")),
        temp_paf=temp("tmp/{sample}_uniform.paf"),
    threads: 1
    resources:
        mem=8,
        hrs=24,
    run:
        with open(input.paf, "r") as infile:
            for i, line in enumerate(infile):
                if i == 0:
                    all_tags = set([x.split(":")[0] for x in line.split("\t")[12:]])
                else:
                    all_tags = all_tags.intersection(
                        set([x.split(":")[0] for x in line.split("\t")[12:]])
                    )

        out_list = []

        with open(input.paf, "r") as infile:
            for line in infile:
                out_list.append(
                    line.split("\t")[0:12]
                    + [x for x in line.split("\t")[12:] if x.split(":")[0] in all_tags]
                )

        with open(output.temp_paf, "w") as outfile:
            for item in out_list:
                outfile.write("\t".join(item))

        df = pd.read_csv(output.temp_paf, sep="\t", low_memory=False, header=None)
        col_out = df.columns.values
        for i, contig in enumerate(df[0].unique()):
            out_num = (i % PARTS) + 1
            df.loc[df[0] == contig][col_out].to_csv(
                f"tmp/{wildcards.sample}.{out_num}-of-{PARTS}.paf",
                sep="\t",
                index=False,
                header=False,
                mode="a+",
            )


rule trim_break_orient_paf:
    input:
        paf="tmp/{sample}.{scatteritem}.paf",
    output:
        contig=temp("tmp/{sample}.{scatteritem}-orient.paf"),
        broken=temp("tmp/{sample}.{scatteritem}-broken.paf"),
    threads: 1
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "rustybam/0.1.27",
    resources:
        mem=24,
        hrs=24,
    shell:
        """
        rustybam orient {input.paf} | rustybam trim-paf | rb filter --paired-len 1000000 > {output.contig}
        rustybam break-paf --max-size {SV_SIZE} {output.contig} > {output.broken}
        """


rule combine_paf:
    input:
        paf=find_contigs,
        flag=rules.split_paf.output.flag,
    output:
        paf="tmp/{sample}-broken.paf",
    threads: 1
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "rustybam/0.1.27",
    resources:
        mem=4,
        hrs=24,
    shell:
        """
        cat {input.paf} > {output.paf} 
        """


rule saff_out:
    input:
        paf=rules.combine_paf.output.paf,
    output:
        saf="results/saffire/{sample}/{sample}.saf",
    threads: 1
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "rustybam/0.1.27",
    resources:
        mem=8,
        hrs=24,
    shell:
        """
        rb stats --paf {input.paf} > {output.saf}
        """
