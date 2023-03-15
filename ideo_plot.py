#!/usr/bin/env python

"""
Demonstrates plotting chromosome ideograms and genes (or any features, really)
using matplotlib.
1) Assumes a file from UCSC's Table Browser from the "cytoBandIdeo" table,
saved as "ideogram.txt". Lines look like this::
    #chrom  chromStart  chromEnd  name    gieStain
    chr1    0           2300000   p36.33  gneg
    chr1    2300000     5300000   p36.32  gpos25
    chr1    5300000     7100000   p36.31  gneg
2) Assumes another file, "ucsc_genes.txt", which is a BED format file
   downloaded from UCSC's Table Browser. This script will work with any
   BED-format file.
"""

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd
import argparse



parser = argparse.ArgumentParser()

parser.add_argument("--ref_name", "-r", type=str, required=False, default='hg38')
parser.add_argument("--asm1", "-a", type=str, required=True)
parser.add_argument("--outfile", "-o", type=str, required=True)
parser.add_argument("--chrom", "-c", type=str, required=False)
parser.add_argument("--sample", "-s", type=str, required=True)
parser.add_argument("--asm2", "-b", type=str, required=False)

args = parser.parse_args()

# args = pd.Series(['chm13', 'results/HG002_rev-1_hifiasm_h1.bed', 'h1.test.pdf'], index=['ref_name', 'asm1', 'outfile'])


# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height, ptype, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('#chrom'):
        # print(chrom)
        xranges = group[['start', 'width']].values
        if ptype != 'ideogram':
            for i, xval in enumerate(xranges):
                if ptype == 'asm1':
                    yield BrokenBarHCollection(
                        [xval], (y_positions[chrom]+(height*(0)**(i%2)), height), facecolors=asm_color_dict[ptype][i%2], **kwargs
                    )
                else:
                    yield BrokenBarHCollection(
                        [xval], (y_positions[chrom]+(height*(0)**(i%2))-(gene_height*2+chrom_height+gene_padding*2), height), facecolors=asm_color_dict[ptype][i%2], **kwargs
                    )
        else:
            yield BrokenBarHCollection(
                    xranges, (y_positions[chrom], height), facecolors=group['colors'], **kwargs
                )
    if del_width:
        del df['width']


asm_color_dict = {'asm1' : {0 : 'tab:blue', 1 : 'tab:orange'}, 'asm2' : {0 : 'tab:green', 1 : 'tab:purple'}}

# Height of each ideogram
chrom_height = 1.25

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 2

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.5


# Spacing between consecutive ideograms
chrom_spacing = (chrom_height+gene_height*4+gene_padding)*1.25
# chrom_spacing = 8

# Width, height (in inches)
figsize = (6, 8)


cyto_dict = {
    'hg38' : '/net/eichler/vol27/projects/hgsvc/nobackups/svpop/data/anno/bands/bands.bed', 
    'chm13' : '/net/eichler/vol26/eee_shared/assemblies/CHM13/T2T/v2.0/anno/cyto.bed'
    }


ideo = pd.read_csv(
    cyto_dict[args.ref_name], 
    sep='\t'
)



# Decide which chromosomes to use
if not args.chrom:
    chromosome_list = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
else:
    chromosome_list = [args.chrom]

ideo = ideo[ideo['#chrom'].apply(lambda x: x in chromosome_list)]

# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0
chrom_ybase = {}
gene_ybase = {}
chrom_centers = {}

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom] = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom] = ybase + chrom_height + gene_padding
    ybase += chrom_height + chrom_spacing


# Add a new column for width
ideo['width'] = ideo['end'] - ideo['start']

# Colors for different chromosome stains
color_lookup = {
    'gneg': (1., 1., 1.),
    'gpos25': (.6, .6, .6),
    'gpos50': (.4, .4, .4),
    'gpos75': (.2, .2, .2),
    'gpos100': (0., 0., 0.),
    'acen': (.8, .4, .4),
    'gvar': (.8, .8, .8),
    'stalk': (.9, .9, .9),
}

# Add a new column for colors
ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])


# Same thing for genes
# df = pd.read_csv(args.asm1, names=['#chrom', 'start', 'end', 'name', 'mapq'])
asm_dict = {}

df = pd.read_csv(args.asm1, names=['#chrom', 'start', 'end', 'name', 'mapq'], sep='\t', header=None)
df = df[df['#chrom'].apply(lambda x: x in chromosome_list)]
df['width'] = df['end'] - df['start']
asm_dict['asm1'] = df.copy()

if args.asm2:
    df = pd.read_csv(args.asm2, names=['#chrom', 'start', 'end', 'name', 'mapq'], sep='\t', header=None)
    df = df[df['#chrom'].apply(lambda x: x in chromosome_list)]
    df['width'] = df['end'] - df['start']
    asm_dict['asm2'] = df.copy()  


fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)


for i, chrom in enumerate(chromosome_list):
    ax.barh(y=chrom_ybase[chrom]+(chrom_height/2), width=max(ideo.loc[ideo['#chrom'] == chrom]['end']), height=chrom_height, edgecolor='black', color="None")

# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, 'ideogram'):
    ax.add_collection(collection)

# ...and the gene data
for asm in asm_dict:
    print(f"plotting {asm}")
    for collection in chromosome_collections(
        asm_dict[asm], gene_ybase, gene_height, asm, alpha=1.0, linewidths=0.5,
    ):
        ax.add_collection(collection)

# Axes tweaking
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list)
ax.set_title(f"{args.sample} vs {args.ref_name.upper()}")
ax.axis('tight')
plt.savefig(args.outfile)
plt.close('all')
