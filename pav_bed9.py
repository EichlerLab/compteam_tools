#!/bin/env python


import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input", "-i", required=True, type=str, help="Input SV bed", nargs='+')
parser.add_argument("--output", "-o", required=True, type=str, help="Output bed9 file")

args = parser.parse_args()



color_dict = {'DEL' : '220,0,0', 'DUP' : '0,0,220', 'INTER' : '0,220,0', 'INV' : '220,140,0', 'TRANSPOSE' : '128,0,128', 'TRANSPOSE_INV' : '128,0,128', 'GAP_NONSYN' : '0,0,0', 'INS' : '0,0,220'}

df = pd.concat(pd.read_csv(x, sep='\t') for x in args.input)

df['COLOR'] = df.apply(lambda row: color_dict[row['SVTYPE']], axis=1)
df['BED_END'] = df.apply(lambda row: row['END']+row['SVLEN'] if row['END']-row['POS'] == 1 else row['END'], axis=1)
df['SCORE'] = '0'
df['STRAND'] = '+'
df['tst'] = 0
df['ten'] = 0
		
df[['#CHROM', 'POS', 'BED_END', 'ID', 'SCORE', 'STRAND', 'tst', 'ten', 'COLOR']].sort_values(['#CHROM', 'POS']).rename(columns={'#CHROM' : '#ct', 'POS' : 'st', 'BED_END' : 'en', 'ID' : 'name', 'SCORE' : 'score', 'STRAND' : 'strand', 'COLOR' : 'color'}).to_csv(args.output, sep='\t', index=False)
