#!/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Generate stacked histograms from an input file and save the plot to an output file.')
parser.add_argument('-i','--input', type=str, help='Input file containing the data', nargs='+')
parser.add_argument('-o', '--output', type=str, help='Output file to save the plot')
args = parser.parse_args()

df = pd.concat([ pd.read_csv(x, sep='\t') for x in args.input ] )

del_small = df.loc[(df['SVTYPE'] == 'DEL') & (df['SVLEN'] <= 500)].copy()
ins_small = df.loc[(df['SVTYPE'] == 'INS') & (df['SVLEN'] <= 500)].copy()
del_large = df.loc[(df['SVTYPE'] == 'DEL') & (df['SVLEN'] < 10000) & (df['SVLEN'] > 500)].copy()
ins_large = df.loc[(df['SVTYPE'] == 'INS') & (df['SVLEN'] < 10000) & (df['SVLEN'] > 500)].copy()
# Define the bins for the histograms
# bins = np.arange(0, 3001, 100)  # Adjust the bin size as needed
bins = 100
# Create a stacked histogram for each plot
plt.figure(figsize=(10, 6))

# Plot for the first set of data
plt.subplot(2, 1, 1)
plt.hist([del_small['SVLEN'], ins_small['SVLEN']], bins=bins, color=['tab:red', 'tab:blue'], alpha=1, label=['DEL', 'INS'], stacked=True)
plt.title('Stacked Histogram of SV Lengths (<= 500bp)')
plt.xlabel('SV Length')
plt.ylabel('Count')
plt.legend()

# Plot for the second set of data
plt.subplot(2, 1, 2)
plt.hist([del_large['SVLEN'], ins_large['SVLEN']], bins=bins, color=['tab:red', 'tab:blue'], alpha=1, label=['DEL', 'INS'], stacked=True)
plt.title('Stacked Histogram of SV Lengths (500bp to 3000bp)')
plt.xlabel('SV Length')
plt.ylabel('Count')
plt.legend()

plt.tight_layout()
plt.savefig(args.output)