#!/bin/env python

import pandas as pd
import os
import numpy as np
import sys
import argparse

sys.path.append('/net/eichler/vol27/projects/structural_variation/nobackups/tools/svpop/202006')

import analib

def find_lim(ylimit):
        # Set scaled y-limit and axis positions
    if ylimit <= 10:
        limits = [0, 5, 10]
        ylimit = 10
    elif ylimit <= 20:
        limits = [0, 10, 20]
        ylimit = 20
    elif ylimit <= 60:
        limits = [0, 30, 60]
        ylimit = 60
    elif ylimit <= 80:
        limits = [0, 40, 80]
        ylimit = 80
    elif ylimit <= 100:
        limits = [0, 50, 100]
        ylimit = 100
    elif ylimit <= 200:
        limits = [0, 100, 200]
        ylimit = 200
    elif ylimit <= 500:
        limits = [0, 250, 500]
        ylimit = 500
    elif ylimit <= 1000:
        limits = [0, 500, 1000]
        ylimit = 1000
    elif ylimit <= 2500:
        limits = [0, 1250, 2500]
        ylimit = 2500
    elif ylimit <= 5000:
        limits = [0, 2500, 5000]
        ylimit = 5000
    elif ylimit <= 7500:
        limits = [0, 3750, 7500]
        ylimit = 7500
    elif ylimit <= 10000:
        limits = [0, 5000, 10000]
        ylimit = 10000
    elif ylimit <= 100000:
        limits = [0, 50000, 100000]
        ylimit = 100000
    elif ylimit <= 150000:
        limits = [0, 75000, 150000]
        ylimit = 150000
    elif ylimit <= 300000:
        limits = [0, 150000, 300000]
        ylimit = 300000
    elif ylimit <= 600000:
        limits = [0, 300000, 600000]
        ylimit = 600000
    elif ylimit <= 1000000:
        limits = [0, 500000, 1000000]
        ylimit = 1000000
    else:
        raise RuntimeError('Ran out of limits: ' + str(ylimit))
    return np.array(limits), ylimit
   

def ideo_cb(df, chrom, ax, fig):
    # Subset to chromosome
    df_a_chrom = df_a.loc[df_a['#CHROM'] == chrom].copy()
    df_b_chrom = df_b.loc[df_b['#CHROM'] == chrom].copy()
    # Bin
    df_a_chrom['BIN_MID'] = ((df_a_chrom['POS'] + df_a_chrom['END']) // 2 // 1e6).astype(np.int32)
    df_b_chrom['BIN_MID'] = ((df_b_chrom['POS'] + df_b_chrom['END']) // 2 // 1e6).astype(np.int32)
    
    bin_max = np.nanmax([
        np.max(df_a_chrom['BIN_MID']) if df_a_chrom.shape[0] > 0 else 0,
        np.max(df_b_chrom['BIN_MID']) if df_b_chrom.shape[0] > 0 else 0
    ])
    if pd.isnull(bin_max):
        bin_max = 0
    x_vals = np.arange(bin_max + 1) * BIN_SIZE  # Inclusive range
    ## Get bar heights ##
    count_hprc_ins = np.zeros(bin_max + 1)
    count_hprc_del = np.zeros(bin_max + 1)
    count_hgsvc_ins = np.zeros(bin_max + 1)
    count_hgsvc_del = np.zeros(bin_max + 1)
    # HPRC INS/DEL
    for val in df_a_chrom.loc[df_a_chrom['SVTYPE'] == 'INS', 'BIN_MID']:
        count_hprc_ins[val] += 1
    for val in df_a_chrom.loc[df_a_chrom['SVTYPE'] == 'DEL', 'BIN_MID']:
        count_hprc_del[val] += 1
    # HGSVC INS/DEL
    for val in df_b_chrom.loc[df_b_chrom['SVTYPE'] == 'INS', 'BIN_MID']:
        count_hgsvc_ins[val] += 1
    for val in df_b_chrom.loc[df_b_chrom['SVTYPE'] == 'DEL', 'BIN_MID']:
        count_hgsvc_del[val] += 1
    # Manually set y-limits (make space for ideo below y=0)
    limits_over, ylim_over = find_lim(np.max((count_hprc_ins + count_hprc_del))*1.05)
    limits_under, ylim_under = find_lim(np.max((count_hgsvc_ins + count_hgsvc_del))*1.05)
    # if ylim_under < ylim_over:
    #     ylim_under_factor = ylim_over/ylim_under
    #     ylim_over_factor = 1
    # else:
    #     ylim_over_factor = ylim_under/ylim_over
    #     ylim_under_factor = 1
    # ideo_space = ylim_over * SPACER_PROP * 2 * ylim_over_factor
    # ylim_high = ylim_over * ylim_over_factor
    # ylim_low = (-ylim_under * ylim_under_factor) - ideo_space
    # count_ideo_space = np.repeat(ideo_space, bin_max + 1)
    # ## Make bar plots ##
    # ax.bar(x_vals, count_hprc_del*ylim_over_factor, width=BIN_SIZE, color='red', label='Deletions')
    # ax.bar(x_vals, count_hprc_ins*ylim_over_factor, width=BIN_SIZE, bottom=count_hprc_del, color='blue', label='Insertions')
    # ax.bar(x_vals, count_hgsvc_del*ylim_under_factor, width=BIN_SIZE, bottom=-(count_hgsvc_del * ylim_under_factor + count_ideo_space), color='red', label=None)
    # ax.bar(x_vals, count_hgsvc_ins*ylim_under_factor, width=BIN_SIZE, bottom=-((count_hgsvc_del + count_hgsvc_ins)*ylim_under_factor + count_ideo_space), color='blue', label=None)
    # ## Set y axis ticks ##
    # ax.set_yticks(
    #     list(np.concatenate([
    #         np.flip(-(limits_under*ylim_under_factor + ideo_space)),
    #         limits_over*ylim_over_factor
    #     ]))
    # )
    # ax.set_yticklabels(
    #     [f'{val:,d}' for val in limits_under[::-1]] + [f'{val:,d}' for val in limits_over]
    # )
    # # Adjust axes
    # ax.set_ylim(ylim_low*ylim_under_factor, ylim_high*ylim_over_factor * (1 + LABEL_SPACE))
    ideo_space = ylim_over * SPACER_PROP * 2 
    ylim_high = ylim_over
    ylim_low = (-ylim_under) - ideo_space
    count_ideo_space = np.repeat(ideo_space, bin_max + 1)
    ## Make bar plots ##
    ax.bar(x_vals, count_hprc_del, width=BIN_SIZE, color='red', label='Deletions')
    ax.bar(x_vals, count_hprc_ins, width=BIN_SIZE, bottom=count_hprc_del, color='blue', label='Insertions')
    ax.bar(x_vals, count_hgsvc_del, width=BIN_SIZE, bottom=-(count_hgsvc_del + count_ideo_space), color='red', label=None)
    ax.bar(x_vals, count_hgsvc_ins, width=BIN_SIZE, bottom=-((count_hgsvc_del + count_hgsvc_ins) + count_ideo_space), color='blue', label=None)
    ## Set y axis ticks ##
    ax.set_yticks(
        list(np.concatenate([
            np.flip(-(limits_under + ideo_space)),
            limits_over
        ]))
    )
    ax.set_yticklabels(
        [f'{val:,d}' for val in limits_under[::-1]] + [f'{val:,d}' for val in limits_over]
    )
    # Adjust axes
    ax.set_ylim(ylim_low, ylim_high * (1 + LABEL_SPACE))




def ideo_mono(df, chrom, ax, fig):
    # Subset to chromosome
    df_all_chrom = df_all.loc[df_all['#CHROM'] == chrom].copy()
    # Bin
    df_all_chrom['BIN_MID'] = ((df_all_chrom['POS'] + df_all_chrom['END']) // 2 // 1e6).astype(np.int32)
    bin_max = np.nanmax([
        np.max(df_all_chrom['BIN_MID']) if df_all_chrom.shape[0] > 0 else 0,
    ])
    if pd.isnull(bin_max):
        bin_max = 0
    x_vals = np.arange(bin_max + 1) * BIN_SIZE  # Inclusive range
    ## Get bar heights ##
    count_hprc_ins = np.zeros(bin_max + 1)
    count_hprc_del = np.zeros(bin_max + 1)
    # HPRC INS/DEL
    for val in df_all_chrom.loc[df_all_chrom['SVTYPE'] == 'INS', 'BIN_MID']:
        count_hprc_ins[val] += 1
    for val in df_all_chrom.loc[df_all_chrom['SVTYPE'] == 'DEL', 'BIN_MID']:
        count_hprc_del[val] += 1
    # HGSVC INS/DEL
    # Manually set y-limits (make space for ideo below y=0)
    # Get limit from data
    limits, ylim = find_lim(np.max(
        [
            np.max((count_hprc_ins + count_hprc_del))
        ]
    ) * 1.05)
    limits = np.array(limits)
    ideo_space = ylim * SPACER_PROP * 2
    ylim_high = ylim
    ylim_low = -ideo_space
    count_ideo_space = np.repeat(ideo_space, bin_max + 1)
    ## Make bar plots ##
    ax.bar(x_vals, count_hprc_del, width=BIN_SIZE, color='red', label='Deletions')
    ax.bar(x_vals, count_hprc_ins, width=BIN_SIZE, bottom=count_hprc_del, color='blue', label='Insertions')
    ## Set y axis ticks ##
    ax.set_yticks(limits)
    ax.set_yticklabels(
        [f'{val:,d}' for val in limits]
    )
    # Adjust axes
    ax.set_ylim(ylim_low, ylim_high * (1 + LABEL_SPACE))



parser = argparse.ArgumentParser() 

parser.add_argument("--a_pattern", "-a", type=str, required=True, help="Bed table with columns ['#CHROM', 'POS', 'END', 'SVTYPE','ID'] where ID == chr-pos-svtype-svlen")
parser.add_argument("--b_pattern", "-b", type=str, required=False, help="Additional bed table if wanting to run a double ideogram or any intersect functions")
parser.add_argument("--intersect", action='store_true', help="If wanting to run a plot of intersecting variants between two samples, requires -i")
parser.add_argument("--only_a", action='store_true', help="If wanting to only plot variants found in A file only")
parser.add_argument("--only_b", action='store_true', help="If wanting to only plot variants found in B file only")
parser.add_argument("--i_file", "-i", type=str, required=False, help="SVPOP intersect file with A,B variants intersected in order")
parser.add_argument("--vartype", "-v", type=str, required=False, default='sv|indel', help="Accepted sv, indel, or sv|indel for both. Default is sv|indel")
parser.add_argument("--svtype", "-s", type=str, required=False, default='ins|del', help="Accepted ins, del, insdel, or ins|del for both. Default is ins|del")
parser.add_argument("--output", "-o", type=str, required=True, help="Output handle for images. Makes both png and pdf")
parser.add_argument("--ref", "-r", type=str, required=False, default='hg38', help="REF to use of ideogram bars: hg38 or chm13_v1.1")

args = parser.parse_args()

patterns = []
intersect = ''

if args.intersect:
    patterns.append(args.a_pattern)
    intersect = 'A,B'   
elif args.only_a:
    patterns.append(args.a_pattern)
    intersect = 'A'
elif args.only_b:
    patterns.append(args.b_pattern)
    intersect = 'B'
elif args.b_pattern:
    patterns.append(args.a_pattern)
    patterns.append(args.b_pattern)
else:
    patterns.append(args.a_pattern)

svtypes = args.svtype.split('|')
vartypes = args.vartype.split('|')

# Paths
# BED_FILT = '/net/eichler/vol27/projects/hprc/nobackups/data_table/qc/hprc_only/hprc_filt/_{svtype}_insdel.bed.gz'

# BED_DROPPED = '/net/eichler/vol27/projects/hprc/nobackups/data_table/qc/hprc_only/hgsvc_filt/{sample}_{svtype}_insdel.bed.gz'

# # BED_ALL = '/net/eichler/vol27/projects/hprc/nobackups/data_table/preqc/hprc_only/tsv/variants_hprc_only_{svtype}_insdel.tsv.gz'
# BED_ALL = '/net/eichler/vol27/projects/marvin_c/nobackups/data_table/cmg_post/tsv/variants_freeze1post_{svtype}_insdel.tsv.gz'

# Fig params
BIN_SIZE = np.int32(1e6)

SPACER_PROP = 0.325  # Shift lower bars by this amount to make space for ideo

LABEL_SPACE = 0.25  # Add this proportion of the y range to the upper limit to make space for the chromosome label


if args.ref == 'chm13_v1.1':
	FAI_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/chm13_v1.1_plus38Y_masked.fasta.fai'
	BAND_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/anno/cyto.bed'
	GAP_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/anno/gap.bed'
	SD_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/anno/sd-max-frac.bed'
	TR_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/anno/trf_regions_200_0.bed'
elif args.ref == 'hg38':
	# Fig tracks
	FAI_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa.fai'
	BAND_FILE_NAME = '/net/eichler/vol27/projects/hgsvc/nobackups/svpop/data/anno/bands/bands.bed'	
	GAP_FILE_NAME = '/net/eichler/vol27/projects/hgsvc/nobackups/svpop/data/anno/gap/gap.bed'	
	SD_FILE_NAME = '/net/eichler/vol27/projects/hgsvc/nobackups/svpop/data/anno/sd/sd-max-frac.bed'
	TR_FILE_NAME = '/net/eichler/vol27/projects/hgsvc/nobackups/svpop/data/anno/trf/trf_regions_200_0.bed'


chroms = ['chr%d' % num for num in range(1,23)]
chroms.append('chrY')
chroms.append('chrX')


if len(patterns) == 1:
    for vartype in vartypes:
        df_all = pd.DataFrame()
        for svtype in svtypes:
            df_all = df_all.append(pd.read_csv(patterns[0].format(svtype=svtype, vartype=vartype), sep='\t'))    
        if intersect == '':
            pass
        else:
            df_int = pd.DataFrame()
            for svtype in ['ins', 'del']:
                df_int = df_int.append(pd.read_csv(args.i_file.format(svtype=svtype, vartype=vartype), sep='\t', header=0))
            if args.intersect:
                df_all = pd.merge(df_all, df_int, left_on='ID', right_on='ID_A')
            else:
                df_all = pd.merge(df_all, df_int, left_on='ID', right_on=f'ID_{intersect}')
            df_all = df_all.loc[df_all['SOURCE_SET'] == intersect]
        ### Figure ###
        # Read ideo bands
        df_band = pd.read_csv(BAND_FILE_NAME, sep='\t')
        df_gap = pd.read_csv(GAP_FILE_NAME, sep='\t')
        df_sd = pd.read_csv(SD_FILE_NAME, sep='\t')
        df_tr = pd.read_csv(TR_FILE_NAME, sep='\t', header=None, names=('#CHROM', 'POS', 'END'))
        # Make figure
        ideo_hist = analib.plot.ideo.ideo_hist(None, FAI_FILE_NAME, df_band, df_gap, df_sd, df_tr, cb_func=ideo_mono)
        # Save
        type_label = "".join(svtypes)
        ideo_hist.fig.savefig(f'{args.output}-{vartype}_{type_label}.png', bbox_inches='tight')
        ideo_hist.fig.savefig(f'{args.output}-{vartype}_{type_label}.pdf', bbox_inches='tight')


if len(patterns) == 2:
    for vartype in vartypes:
        df_a = pd.DataFrame()
        df_b = pd.DataFrame()
        for svtype in svtypes:
            df_a = df_a.append(pd.read_csv(patterns[0].format(svtype=svtype, vartype=vartype), sep='\t', usecols=['#CHROM', 'POS', 'END', 'SVTYPE', 'ID']).dropna())
            df_b = df_b.append(pd.read_csv(patterns[1].format(svtype=svtype, vartype=vartype), sep='\t',usecols=['#CHROM', 'POS', 'END', 'SVTYPE', 'ID']).dropna())
        ### Figure ###
        # Read ideo bands
        df_band = pd.read_csv(BAND_FILE_NAME, sep='\t')
        df_gap = pd.read_csv(GAP_FILE_NAME, sep='\t')
        df_sd = pd.read_csv(SD_FILE_NAME, sep='\t')
        df_tr = pd.read_csv(TR_FILE_NAME, sep='\t', header=None, names=('#CHROM', 'POS', 'END'))
        # Make figure
        ideo_hist = analib.plot.ideo.ideo_hist(None, FAI_FILE_NAME, df_band, df_gap, df_sd, df_tr, cb_func=ideo_cb)
        # Save
        type_label = "".join(svtypes)
        ideo_hist.fig.savefig(f'{args.output}-{vartype}_{type_label}.png', bbox_inches='tight')
        ideo_hist.fig.savefig(f'{args.output}-{vartype}_{type_label}.pdf', bbox_inches='tight')

