#!/bin/env python

import pandas as pd
import argparse
import os, sys
from datetime import datetime
import hashlib
import logging

def _get_checksum(file_name):
    """
    Get the MD5 checksum of a file.
    :param file_name: File to check.
    :return: MD5 checksum.
    """
    # Code from Stack Overflow:
    # http://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file                                                            
    hash_md5 = hashlib.md5()
    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


logging.basicConfig(
    format="%(levelname)s (%(asctime)s): %(message)s (Line: %(lineno)d [%(filename)s])",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.INFO,
)


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter) 

parser.add_argument('--input', '-i', type=str, required=True, help='Input copy record')
parser.add_argument('--sample_new', '-s', type=str, required=True, help='Correct sample')
parser.add_argument('--cohort', '-c', type=str, required=False, help='Cohort of correct sample')
parser.add_argument('--output', '-o', type=str, required=False, help='Output file with files to remove')

args = parser.parse_args()


if not args.output:
    if "PacBio_HiFi" in args.input:
        output_file = "-".join([args.sample_new, "REMOVE"] + args.input.split("/")[-2:])
else:
    output_file = args.output

output_file = output_file.replace(".gz", "")

df = pd.read_csv(args.input, sep="\t")

new_files = pd.DataFrame()

for index in df.index:
    # added notna() condition to prevent the error from .endswith with null value in DEST_PATH
    if pd.notna(df.iloc[index]['DEST_PATH']) and df.iloc[index]['DEST_PATH'].endswith("bam") and os.path.isfile(df.iloc[index]['DEST_PATH'].replace("bam", "fastq.gz.fai")):
        logging.info("Finding fastq files associated with "+df.iloc[index]['DEST_PATH'])
        for file_end in ["fastq.gz", "fastq.gz.fai", "fastq.gz.gzi"]:
            new_file = df.iloc[index]['DEST_PATH'].replace("bam", file_end)
            if new_file in df['DEST_PATH'].values:
                continue
            file_stat = os.stat(new_file)
            md5 = _get_checksum(new_file)
            file_dict = {
                        "SAMPLE" : ["NA" ],
                        "SEQ_TYPE" : [df.iloc[index]["SEQ_TYPE"]],
                        "RUN_ID": [df.iloc[index]["RUN_ID"]],
                        "CELL" : [df.iloc[index]["CELL"]],
                        "CCS_JOB_ID" : ["MANUAL"],
                        "BARCODE" : [df.iloc[index]["BARCODE"]],
                        "SOURCE_PATH" : ["NA"],
                        "DEST_PATH" : [new_file],
                        "SIZE" : [file_stat.st_size],
                        "MOD_TIME" : [datetime.fromtimestamp(file_stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S.%f")],
                        "STATUS" : ["Copied"],
                        "MD5" : [md5],
                        }
            new_files = pd.concat([new_files, pd.DataFrame.from_dict(file_dict)])

df  = pd.concat([df, new_files]).reset_index(drop=True).fillna("NA")

if args.cohort:
    new_prefix = os.path.join(f"../{args.cohort}", args.sample_new)
else:
    new_prefix = args.sample_new

df['COPY_PATH'] = df["DEST_PATH"].apply(lambda val : os.path.join(new_prefix, "/".join(val.split("/")[1:])))
df['NEW_PATH'] = df["DEST_PATH"].apply(lambda val : os.path.join(args.sample_new, "/".join(val.split("/")[1:])))

logging.info("Linking files")
migrated_files = [args.input]
for index in df.index:
    if os.path.exists(df.iloc[index]['DEST_PATH']):
        os.makedirs(os.path.dirname(df.iloc[index]['COPY_PATH']), exist_ok=True)
        # added to aviod file exist Error
        if os.path.exists(df.iloc[index]['COPY_PATH']):
            os.remove(df.iloc[index]['COPY_PATH'])
        os.link(df.iloc[index]['DEST_PATH'], df.iloc[index]['COPY_PATH'])
        migrated_files.append(df.iloc[index]['DEST_PATH'])

run_id = df.iloc[0]['RUN_ID']
try:
    bc = df.iloc[0]['BC']
except:
    bc = df.iloc[0]['BARCODE']
cell = df.iloc[0]['CELL']
old_sample = args.input.split("/")[0]


for tsv in ["cell_table.tsv", "cell_table_revio.tsv"]:
    cell_df = pd.read_csv(tsv, sep="\t")
    cell_df = cell_df.loc[~((cell_df["BARCODE"] == bc) & (cell_df["SAMPLE"] == old_sample) & (cell_df["RUN_ID"] == "RUN_ID") & (cell_df["CELL"] == cell))]
    cell_df.to_csv(tsv, sep="\t", index=False)


df['DEST_PATH'] = df['NEW_PATH']


out_cols = [
"SAMPLE",
"SEQ_TYPE",
"RUN_ID",
"CELL",
"CCS_JOB_ID",
"BARCODE",
"SOURCE_PATH",
"DEST_PATH",
"SIZE",
"MOD_TIME",
"STATUS",
"MD5"
]


run_id = os.path.basename(args.input).replace(".tab.gz", "")
old_sample = args.input.split("/")[0]

copy_rename = os.path.join(new_prefix, "/".join(args.input.split("/")[1:]))

os.makedirs(os.path.dirname(copy_rename), exist_ok=True)

df[out_cols].to_csv(copy_rename, sep="\t", index=False)

with open(output_file, "w") as outfile:
    outfile.write("\n".join(migrated_files)+"\n")

logging.info(f"Output written to {output_file}")
