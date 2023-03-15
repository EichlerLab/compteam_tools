#!/usr/bin/env python3
"""
Filter PAV calls (only supports SV-DEL, SV-INS, SV-INV)
Usage: ./pav_filter_more.py --query path/to/asm/bed/sv_del.bed.gz --sample 14455_p1
"""

from pybedtools import BedTool
import pandas as pd
import os
import argparse
import sys

ANNO_DIR = '/net/eichler/vol27/projects/autism_genome_assembly/nobackups/anno'
PAV_DIR = '/net/eichler/vol27/projects/autism_genome_assembly/nobackups/analysis/pav-trios'
ASSEMBLY_EVAL_DIR = '/net/eichler/vol27/projects/autism_genome_assembly/nobackups/analysis/assembly_eval'

lc_bed = os.path.join(ANNO_DIR, "lc.bed.gz")
pli_bed = os.path.join(ANNO_DIR, "pliByGene-052022.bed.gz")
exon_bed = os.path.join(ANNO_DIR, "refseq_curated-grch38_exons.bed.gz")

single_header = ['contig', 'start', 'end', 'name', 'gene', 'leouf', 'exon_overlaps', 'genotype', 'ae_flagged']
trio_header = ['seen_in_parents']
quad_header = ['seen_in_sibling'] + trio_header


def filter_sizes(query_bed, size=30000):
    # NOTE: pav SVLEN falls in the strand column and SVTYPE falls in score column
    size_filtered = BedTool(query_bed).filter(lambda x: int(x.strand) >= size)

    def add_len(feature):
        if 'INS' in feature.score:
            feature.end = int(feature.end) + int(feature.strand)
            return feature
        else:
            return feature

    return size_filtered.each(add_len)


def filter_more(query_bed):
    # remove low-complexity regions
    a = query_bed.intersect(lc_bed, v=True)

    # overlap with leouf/pli scores
    b = a.intersect(pli_bed, wa=True, wb=True)

    # overlap with exon
    overlaps = b.intersect(exon_bed, wa=True, wb=True)
    df = overlaps.to_dataframe(usecols=[0, 1, 2, 3, 7, 10, 25, 28, 32],
                               names=['contig', 'start', 'end', 'name', 'genotype', 'query_space', 'gene', 'leouf',
                                      'exon'])
    df['exon_overlaps'] = df.groupby(['name', 'gene'])['contig'].transform('size')
    df.drop_duplicates(subset='gene', inplace=True)

    df.sort_values(by=['leouf'], ascending=True, inplace=True)
    return df


def seen_in_parents(sample, regions, variant_type='sv_ins') -> dict:
    family = sample.split('_')[0]
    parents = ['mo', 'fa']
    seen = {}
    for member in parents:
        family_id = f'{family}_{member}'
        bed_path = os.path.join(PAV_DIR, 'results', family_id, f'bed/{variant_type}.bed.gz')
        bed = BedTool(regions).intersect(bed_path)
        collect_intervals = []
        for interval in bed:
            collect_intervals.append(concat_fields(interval))
        seen[member] = collect_intervals
    return seen


def seen_in_sibling(sample, regions, variant_type='sv_ins') -> list:
    bed_path = os.path.join(PAV_DIR, 'results', sample, f'bed/{variant_type}.bed.gz')
    seen = []
    bed = BedTool(regions).intersect(bed_path)
    for interval in bed:
        seen.append(concat_fields(interval))
    return seen


def concat_fields(interval) -> str:
    """
    Concatenate fields to use as index in a pandas dataframe.
    :param interval: BedTools interval object
    :return: chrom:start-end
    """
    return '{0}:{1}-{2}'.format(interval['chrom'], interval['start'], interval['end'])


def flagged_regions(sample, regions, by='assembly_eval') -> list:
    """
    Check if regions are flagged by assembly_eval.
    :param sample: e.g. 14455_p1
    :param regions: a list of strings e.g. ['hg1:1-10','hg1:10-11']
    :param by:
    :return: a list of flagged regions
    """
    # only take one haplotype if it's found in both - FOR NOW
    new_regions = [x.split(';')[0] if ';' in x else x for x in regions]
    df = pd.DataFrame(new_regions, columns=['query_space'])
    df[['contig', 'st', 'end']] = df['query_space'].str.split(r'\W', expand=True)
    df.set_index('query_space', inplace=True)
    query_bed = BedTool.from_dataframe(df)
    seen = []
    qc_bed_path = os.path.join(ASSEMBLY_EVAL_DIR, sample, 'assembly_eval/HiFi/minimap2/final_qc.bed')
    if os.path.isfile(qc_bed_path):
        bed = query_bed.intersect(qc_bed_path)
        for interval in bed:
            seen.append(concat_fields(interval))
    else:
        raise FileNotFoundError(qc_bed_path)
    return seen


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser()

    parser.add_argument('--filter', type=int, default=30000,
                        help='Size filter')
    parser.add_argument('--trio', default=False, required=False,
                        action='store_true', help='If the ASM/sample is part of a trio.')
    parser.add_argument('--quad', type=str, default=False, required=False,
                        help='Supply 1 sibling id e.g. 14455_s1')
    parser.add_argument('--anno_dir', required=False,
                        help='Directory holding the database annotation files.')
    parser.add_argument('--output', '-o', type=str, required=False, help='Output file to write to.',
                        default='/dev/stdout')

    required_named = parser.add_argument_group('required arguments')
    required_named.add_argument('--query', type=str, required=True,
                        help='PAV bed file')
    required_named.add_argument('--sample', type=str, required=True,
                        help='e.g. 14455_p1')
    return parser


def main():
    """Run script to filter PAV bed files"""
    parser = get_parser()
    args = parser.parse_args()

    # filter by size
    size_filtered = filter_sizes(args.query, args.filter)

    # annotate with leouf and exons
    target_df = filter_more(size_filtered)

    # annotate presence in assembly_eval pipeline
    ae_flagged = flagged_regions(regions=target_df['query_space'].to_list(), sample=args.sample)
    if not ae_flagged:
        target_df['ae_flagged'] = ['no'] * len(target_df)
    else:
        target_df['ae_flagged'] = ['yes' if x in ae_flagged else 'no' for x in target_df['query_space'].values]

    # extract sv_type from path name; e.g. sv_ins
    variant = os.path.basename(args.query).split('.')[0]
    
    if args.trio:
        temp_bed = BedTool.from_dataframe(target_df)
        target_df.index = [concat_fields(interval) for interval in temp_bed]
        parents_dict = seen_in_parents(regions=temp_bed, sample=args.sample, variant_type=variant)
        del temp_bed
        target_df['seen_in_mo'] = ['mo' if name in parents_dict['mo'] else '' for name in target_df.index]
        target_df['seen_in_fa'] = ['fa' if name in parents_dict['fa'] else '' for name in target_df.index]
        target_df['seen_in_parents'] = target_df["seen_in_mo"] + ',' + target_df["seen_in_fa"]
        print(target_df.to_csv(sep='\t', index=False, columns=single_header + trio_header))
    elif args.quad:
        temp_bed = BedTool.from_dataframe(target_df)
        target_df.index = [concat_fields(interval) for interval in temp_bed]
        sibling = seen_in_sibling(regions=temp_bed, sample=args.quad, variant_type=variant)
        parents_dict = seen_in_parents(regions=temp_bed, sample=args.sample, variant_type=variant)
        del temp_bed
        target_df['seen_in_mo'] = ['mo' if name in parents_dict['mo'] else '' for name in target_df.index]
        target_df['seen_in_fa'] = ['fa' if name in parents_dict['fa'] else '' for name in target_df.index]
        target_df['seen_in_parents'] = target_df["seen_in_mo"] + ',' + target_df["seen_in_fa"]
        target_df['seen_in_sibling'] = ['yes' if name in sibling else 'no' for name in target_df.index]
        print(target_df.to_csv(sep='\t', index=False, columns=single_header + quad_header))
    else:
        print(target_df.to_csv(sep='\t', index=False, columns=single_header))


if __name__ == "__main__":
    sys.exit(main())
