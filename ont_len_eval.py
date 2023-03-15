#!/bin/env python3

import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser() 

parser.add_argument("--input", "-i", type=str, required=False, help="Input file name")
parser.add_argument("--fofn", "-f", type=str, required=False, help="Input fofn file name")
parser.add_argument("--output", "-o", type=str, required=False, default="/dev/stdout", help="Output file name")
parser.add_argument("--denom", "-d", type=int, required=False, default=1000000000, help="Denominator to divide raw bases by (default: 1000000000}")
parser.add_argument("--cutoffs", "-c", type=str, required=False, default="80,100", help="Threshold for summing in kb (default: 80,100)")

args = parser.parse_args()


if args.input:
	df = pd.read_csv(args.input, sep='\t')
	len_df = df.loc[df['passes_filtering'] == True]
	with open(args.output, 'w') as outFile:
		len_columns = '\t'.join([f'cov >{x} kbp' for x in args.cutoffs.split(',')])
		outFile.write(f'sample\ttotal cov\t{len_columns}\n')
		len_values = '\t'.join([str(np.sum(len_df.loc[len_df['sequence_length_template'] > int(x)*1000]['sequence_length_template'])/float(args.denom)) for x in args.cutoffs.split(',')])
		sample = len_df.iloc[0]['sample_id']
		total_cove = np.sum(len_df['sequence_length_template'])/float(args.denom)
		outFile.write(f'{sample}\t{total_cove}\t{len_values}\n')
elif args.fofn:
	with open(args.output, 'w') as outFile:
		len_columns = '\t'.join([f'cov >{x} kbp' for x in args.cutoffs.split(',')])
		outFile.write(f'sample\ttotal cov\t{len_columns}\n')
		with open(args.fofn, 'r') as inFile:
			for line in inFile:
				df = pd.read_csv(line.rstrip(), sep='\t')
				len_df = df.loc[df['passes_filtering'] == True]
				len_values = '\t'.join([str(np.sum(len_df.loc[len_df['sequence_length_template'] > int(x)*1000]['sequence_length_template'])/float(args.denom)) for x in args.cutoffs.split(',')])
				sample = len_df.iloc[0]['sample_id']
				total_cove = np.sum(len_df['sequence_length_template'])/float(args.denom)
				outFile.write(f'{sample}\t{total_cove}\t{len_values}\n')
else:
	print('Either -f or -i must be declared')

