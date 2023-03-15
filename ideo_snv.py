#!/bin/env python

import pandas as pd
import os
import numpy as np
import sys
import argparse

sys.path.append('/net/eichler/vol27/projects/structural_variation/nobackups/tools/svpop/202006')

import analib

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
    count_a_snv = np.zeros(bin_max + 1)
    count_b_snv = np.zeros(bin_max + 1)
    # HPRC INS/DEL
    for val in df_a_chrom.loc[df_a_chrom['SVTYPE'] == 'SNV', 'BIN_MID']:
        count_a_snv[val] += 1
    # HGSVC INS/DEL
    for val in df_b_chrom.loc[df_b_chrom['SVTYPE'] == 'SNV', 'BIN_MID']:
        count_b_snv[val] += 1
    # Manually set y-limits (make space for ideo below y=0)
    # Get limit from data
    ylim = np.max(
        [
            np.max(count_a_snv),
            np.max(count_b_snv)
        ]
    ) * 1.05
    # Set scaled y-limit and axis positions
    if ylim <= 10:
        limits = [0, 5, 10]
        ylim = 10
    elif ylim <= 20:
        limits = [0, 10, 20]
        ylim = 20
    elif ylim <= 60:
        limits = [0, 30, 60]
        ylim = 60
    elif ylim <= 80:
        limits = [0, 40, 80]
        ylim = 80
    elif ylim <= 100:
        limits = [0, 50, 100]
        ylim = 100
    elif ylim <= 200:
        limits = [0, 100, 200]
        ylim = 200
    elif ylim <= 500:
        limits = [0, 250, 500]
        ylim = 500
    elif ylim <= 1000:
        limits = [0, 500, 1000]
        ylim = 1000
    elif ylim <= 2500:
        limits = [0, 1250, 2500]
        ylim = 2500
    elif ylim <= 5000:
        limits = [0, 2500, 5000]
        ylim = 5000
    elif ylim <= 7500:
        limits = [0, 3750, 7500]
        ylim = 7500
    elif ylim <= 10000:
        limits = [0, 5000, 10000]
        ylim = 10000
    elif ylim <= 100000:
        limits = [0, 50000, 100000]
        ylim = 100000
    elif ylim <= 150000:
        limits = [0, 75000, 150000]
        ylim = 150000
    elif ylim <= 300000:
        limits = [0, 150000, 300000]
        ylim = 300000
    elif ylim <= 600000:
        limits = [0, 300000, 600000]
        ylim = 600000
    elif ylim <= 1000000:
        limits = [0, 500000, 1000000]
        ylim = 1000000
    else:
        raise RuntimeError('Ran out of limits: ' + str(ylim))
    limits = np.array(limits)
    ideo_space = ylim * SPACER_PROP * 2
    ylim_high = ylim
    ylim_low = -ylim - ideo_space
    count_ideo_space = np.repeat(ideo_space, bin_max + 1)
    ## Make bar plots ##
    ax.bar(x_vals, count_a_snv, width=BIN_SIZE, color=args.color, label='Deletions')
    ax.bar(x_vals, count_b_snv, width=BIN_SIZE, bottom=-(count_a_snv + count_ideo_space), color=args.color, label=None)
    ## Set y axis ticks ##
    ax.set_yticks(
        list(np.concatenate([
            np.flip(-(limits + ideo_space)),
            limits
        ]))
    )
    ax.set_yticklabels(
        [f'{val:,d}' for val in limits[::-1]] + [f'{val:,d}' for val in limits]
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
    count_snv = np.zeros(bin_max + 1)
    # HPRC INS/DEL
    for val in df_all_chrom.loc[df_all_chrom['SVTYPE'] == 'SNV', 'BIN_MID']:
        count_snv[val] += 1
    # HGSVC INS/DEL
    # Manually set y-limits (make space for ideo below y=0)
    # Get limit from data
    ylim = np.max(
        count_snv
    ) * 1.05
    # Set scaled y-limit and axis positions
    if ylim <= 10:
        limits = [0, 5, 10]
        ylim = 10
    elif ylim <= 20:
        limits = [0, 10, 20]
        ylim = 20
    elif ylim <= 60:
        limits = [0, 30, 60]
        ylim = 60
    elif ylim <= 80:
        limits = [0, 40, 80]
        ylim = 80
    elif ylim <= 100:
        limits = [0, 50, 100]
        ylim = 100
    elif ylim <= 200:
        limits = [0, 100, 200]
        ylim = 200
    elif ylim <= 500:
        limits = [0, 250, 500]
        ylim = 500
    elif ylim <= 1000:
        limits = [0, 500, 1000]
        ylim = 1000
    elif ylim <= 2500:
        limits = [0, 1250, 2500]
        ylim = 2500
    elif ylim <= 5000:
        limits = [0, 2500, 5000]
        ylim = 5000
    elif ylim <= 7500:
        limits = [0, 3750, 7500]
        ylim = 7500
    elif ylim <= 10000:
        limits = [0, 5000, 10000]
        ylim = 10000
    elif ylim <= 100000:
        limits = [0, 50000, 100000]
        ylim = 100000
    elif ylim <= 150000:
        limits = [0, 75000, 150000]
        ylim = 150000
    elif ylim <= 300000:
        limits = [0, 150000, 300000]
        ylim = 300000
    elif ylim <= 600000:
        limits = [0, 300000, 600000]
        ylim = 600000
    elif ylim <= 1000000:
        limits = [0, 500000, 1000000]
        ylim = 1000000
    else:
        raise RuntimeError('Ran out of limits: ' + str(ylim))
    limits = np.array(limits)
    ideo_space = ylim * SPACER_PROP * 2
    ylim_high = ylim
    ylim_low = -ideo_space
    count_ideo_space = np.repeat(ideo_space, bin_max + 1)
    ## Make bar plots ##
    ax.bar(x_vals, count_snv, width=BIN_SIZE, color=args.color, label='SNVs')
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
parser.add_argument("--output", "-o", type=str, required=True, help="Output handle for images. Makes both png and pdf")
parser.add_argument("--color", "-c", type=str, required=False, default='tab:green', help="Color to plot SNVs with")

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


# Paths
# BED_FILT = '/net/eichler/vol27/projects/hprc/nobackups/data_table/qc/hprc_only/hprc_filt/_{svtype}_insdel.bed.gz'

# BED_DROPPED = '/net/eichler/vol27/projects/hprc/nobackups/data_table/qc/hprc_only/hgsvc_filt/{sample}_{svtype}_insdel.bed.gz'

# # BED_ALL = '/net/eichler/vol27/projects/hprc/nobackups/data_table/preqc/hprc_only/tsv/variants_hprc_only_{svtype}_insdel.tsv.gz'
# BED_ALL = '/net/eichler/vol27/projects/marvin_c/nobackups/data_table/cmg_post/tsv/variants_freeze1post_{svtype}_insdel.tsv.gz'

# Fig params
BIN_SIZE = np.int32(1e6)

SPACER_PROP = 0.325  # Shift lower bars by this amount to make space for ideo

LABEL_SPACE = 0.25  # Add this proportion of the y range to the upper limit to make space for the chromosome label

# Fig tracks
FAI_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/chm13_v1.1_plus38Y_masked.fasta.fai'

BAND_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/anno/cyto.bed'
GAP_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/anno/gap.bed'
SD_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/anno/sd-max-frac.bed'
TR_FILE_NAME = '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v1.1/anno/trf_regions_200_0.bed'


chroms = ['chr%d' % num for num in range(1,23)]
chroms.append('chrY')
chroms.append('chrX')


if len(patterns) == 1:
    df_all = pd.read_csv(patterns[0], sep='\t')
    if intersect == '':
        df_band = pd.read_csv(BAND_FILE_NAME, sep='\t')
        df_gap = pd.read_csv(GAP_FILE_NAME, sep='\t')
        df_sd = pd.read_csv(SD_FILE_NAME, sep='\t')
        df_tr = pd.read_csv(TR_FILE_NAME, sep='\t', header=None, names=('#CHROM', 'POS', 'END'))
        # Make figure
        ideo_hist = analib.plot.ideo.ideo_hist(None, FAI_FILE_NAME, df_band, df_gap, df_sd, df_tr, cb_func=ideo_mono)
        # Save
        ideo_hist.fig.savefig(f'{args.output}-snv_snv.png', bbox_inches='tight')
        ideo_hist.fig.savefig(f'{args.output}-snv_snv.pdf', bbox_inches='tight')
    else:
        df_int = pd.read_csv(args.i_file, sep='\t', header=0)
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
        ideo_hist.fig.savefig(f'{args.output}-snv_snv.png', bbox_inches='tight')
        ideo_hist.fig.savefig(f'{args.output}-snv_snv.pdf', bbox_inches='tight')


if len(patterns) == 2:
    df_a = pd.read_csv(patterns[0], sep='\t', usecols=['#CHROM', 'POS', 'END', 'SVTYPE', 'ID']).dropna()
    df_b = pd.read_csv(patterns[1], sep='\t', usecols=['#CHROM', 'POS', 'END', 'SVTYPE', 'ID']).dropna()
    ### Figure ###
    # Read ideo bands
    df_band = pd.read_csv(BAND_FILE_NAME, sep='\t')
    df_gap = pd.read_csv(GAP_FILE_NAME, sep='\t')
    df_sd = pd.read_csv(SD_FILE_NAME, sep='\t')
    df_tr = pd.read_csv(TR_FILE_NAME, sep='\t', header=None, names=('#CHROM', 'POS', 'END'))
    # Make figure
    ideo_hist = analib.plot.ideo.ideo_hist(None, FAI_FILE_NAME, df_band, df_gap, df_sd, df_tr, cb_func=ideo_cb)
    # Save
    ideo_hist.fig.savefig(f'{args.output}-snv_snv.png', bbox_inches='tight')
    ideo_hist.fig.savefig(f'{args.output}-snv_snv.pdf', bbox_inches='tight')

