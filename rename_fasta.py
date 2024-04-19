#!/bin/env python3

import argparse
from Bio import SeqIO

def rename_fasta_headers(input_file, output_file, sample_name):
    """
    Rename the headers of a FASTA file sequentially.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file with renamed headers.
    """
    # Open input and output files
    seq_dict = {}
    with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
        # Iterate through each record in the input FASTA file
        for i, record in enumerate(SeqIO.parse(f_in, "fasta"), 1):
            # Generate a new header using the record's index
            new_header = f"{sample_name}#CTG{i:04}"
            seq_dict[record.name] = new_header
            # Modify the header of the record
            record.description = record.name = new_header
            # Write the modified record to the output FASTA file
            SeqIO.write(record, f_out, "fasta")
    return seq_dict

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Rename FASTA headers.")
    parser.add_argument("input_file", help="Path to the input FASTA file")
    parser.add_argument("output_file", help="Path to the output FASTA file")
    parser.add_argument("sample_name", help="Sample name to be used in fasta headers")
    parser.add_argument("rename_tsv", help="Path to write out the rename TSV")
    args = parser.parse_args()

    # Rename FASTA headers
    seq_dict = rename_fasta_headers(args.input_file, args.output_file, args.sample_name)
    with open(args.rename_tsv, "w") as outfile:
        for orig_name in seq_dict:
            outfile.write("\n".join("\t".join([orig_name, seq_conv_all[orig_name]])))



if __name__ == "__main__":
    main()
