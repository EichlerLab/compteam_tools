#!/bin/env python3

import argparse
from Bio import SeqIO
import gzip

def opener(filename) -> object:
    # Reads in first two bytes of the input file and looks for gzip bytes
    # Returns gzip.open if found, else open
    f = open(filename, "rb")
    if f.read(2) == b"\x1f\x8b":
        f.close()
        return gzip.open
    else:
        f.close()
        return open


def find_n_ranges(fasta_file, output_bed):
    # Open the output BED file for writing
    with open(output_bed, 'w') as bed_file:
        open_func = opener(fasta_file)
        with open_func(fasta_file, "rt") as file_in:
            # Parse the FASTA file and iterate over each sequence
            for record in SeqIO.parse(file_in, 'fasta'):
                seq_id = record.id
                sequence = str(record.seq).upper()
                
                n_ranges = []
                start = None

                # Iterate through the sequence, identifying 'N' ranges
                for i, base in enumerate(sequence):
                    if base == args.base:
                        if start is None:
                            start = i
                    else:
                        if start is not None:
                            n_ranges.append((start, i - 1))
                            start = None

                # Check if there is an 'N' range at the end of the sequence
                if start is not None:
                    n_ranges.append((start, len(sequence) - 1))

                # Write the 'N' ranges to the BED file
                for start, end in n_ranges:
                    bed_file.write(f"{seq_id}\t{start}\t{end}\n")

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate a BED file with ranges of the letter N from a FASTA file.')
parser.add_argument('input_fasta', help='Path to the input FASTA file')
parser.add_argument('output_bed', help='Path to the output BED file')
parser.add_argument('--base', '-b', help='Base to count (default: N)', default='N')

# Retrieve the input and output file paths from command line arguments
args = parser.parse_args()
input_fasta = args.input_fasta
output_bed = args.output_bed

# Call the function to generate the BED file with 'N' ranges
find_n_ranges(input_fasta, output_bed)
