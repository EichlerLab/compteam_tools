#!/bin/env python

import pandas as pd
import argparse
import os, sys
from datetime import datetime
import hashlib
import logging
import glob
from pathlib import Path

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
    output_file =  "-".join([args.sample_new, "REMOVE"] + [args.input.split("/")[-1]])
else:
    output_file = args.output

output_file = output_file.replace(".gz", "")

df = pd.read_csv(args.input, sep="\t")


new_files = pd.DataFrame()

logging.info("Looking for fast5 and pod5, this might take a while")
for index in df.index:
    if df.iloc[index]['DEST_PATH'].endswith("fast5"):
        search_base = "/".join(df.iloc[index]['DEST_PATH'].split("/")[0:-2])
        search_folder = df.iloc[index]['DEST_PATH'].split("/")[-2]
        search_folder_pod = df.iloc[index]['DEST_PATH'].split("/")[-2].replace('fast5', 'pod5')
        search_file = df.iloc[index]['DEST_PATH'].split("/")[-1]
        search_file_pod = df.iloc[index]['DEST_PATH'].split("/")[-1].replace('fast5', 'pod5')
        for folder_check in [search_folder, search_folder_pod]:
            for file_check in [search_file, search_file_pod]:
                if os.path.isfile(os.path.join(search_base, folder_check, file_check)):
                    new_file = os.path.join(search_base, folder_check, file_check)
                    if new_file in df['DEST_PATH'].values:
                        continue
                    logging.debug(f"Getting info for NEW_FILE: {new_file}")
                    file_stat = os.stat(new_file)
                    md5 = _get_checksum(new_file)
                    file_dict = {
                                "SAMPLE" : ["NA" ],
                                "SEQ_TYPE" : [df.iloc[index]["SEQ_TYPE"]],
                                "RUN_ID": [df.iloc[index]["RUN_ID"]],
                                "SOURCE_PATH" : ["NA"],
                                "DEST_PATH" : [new_file],
                                "SIZE" : [file_stat.st_size],
                                "MOD_TIME" : [datetime.fromtimestamp(file_stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S.%f")],
                                "STATUS" : ["Copied"],
                                "MD5" : [md5],
                                }
                    new_files = pd.concat([new_files, pd.DataFrame.from_dict(file_dict)])

fastq_files = pd.DataFrame()

if os.path.exists(os.path.join("/".join(args.input.split("/")[0:4]), "fastq", os.path.basename(args.input).replace(".tab.gz", ""))):
    fastq_dir = os.path.join("/".join(args.input.split("/")[0:4]), "fastq", os.path.basename(args.input).replace(".tab.gz", ""))
    for path in Path(fastq_dir).rglob("*"):
        if os.path.isfile(os.path.join(path.parent, path.name)):
            fastq_files = pd.concat([fastq_files, pd.DataFrame.from_dict({"DEST_PATH" : [os.path.join(path.parent, path.name)]})])


fastq_files = fastq_files.reset_index(drop=True)

df  = pd.concat([df, new_files]).reset_index(drop=True).fillna("NA")

if args.cohort:
    new_prefix = os.path.join(f"../{args.cohort}", args.sample_new)
else:
    new_prefix = args.sample_new

old_sample = args.input.split("/")[0]


df['COPY_PATH'] = df["DEST_PATH"].apply(lambda val : os.path.join(new_prefix, "/".join(val.split("/")[1:])))
df['NEW_PATH'] = df["DEST_PATH"].apply(lambda val : os.path.join(args.sample_new, "/".join(val.split("/")[1:])))

if len(fastq_files) > 0:
    fastq_files['COPY_PATH'] = fastq_files["DEST_PATH"].apply(lambda val : os.path.join(new_prefix, "/".join(val.split("/")[1:])).replace(old_sample, args.sample_new))
    fastq_files['NEW_PATH'] = fastq_files["DEST_PATH"].apply(lambda val : os.path.join(args.sample_new, "/".join(val.split("/")[1:])).replace(old_sample, args.sample_new))

logging.info("Linking files")
migrated_files = [args.input]
for index in df.index:
    if os.path.exists(df.iloc[index]['DEST_PATH']):
        os.makedirs(os.path.dirname(df.iloc[index]['COPY_PATH']), exist_ok=True)
        os.link(df.iloc[index]['DEST_PATH'], df.iloc[index]['COPY_PATH'])
        migrated_files.append(df.iloc[index]['DEST_PATH'])


for index in fastq_files.index:
    if os.path.exists(fastq_files.iloc[index]['DEST_PATH']):
        os.makedirs(os.path.dirname(fastq_files.iloc[index]['COPY_PATH']), exist_ok=True)
        os.link(fastq_files.iloc[index]['DEST_PATH'], fastq_files.iloc[index]['COPY_PATH'])
        migrated_files.append(fastq_files.iloc[index]['DEST_PATH'])

df['SAMPLE'] = args.sample_new
df['DEST_PATH'] = df['NEW_PATH']


out_cols = [
"SAMPLE",
"SEQ_TYPE",
"RUN_ID",
"SOURCE_PATH",
"DEST_PATH",
"SIZE",
"MOD_TIME",
"STATUS",
"MD5"
]

run_id = os.path.basename(args.input).replace(".tab.gz", "")
old_sample = args.input.split("/")[0]


logging.info("Rewriting fast5 and basecall tsv files")
fast5_df = pd.read_csv("ont_fast5_table.tsv", sep='\t')
fast5_df = fast5_df.loc[~((fast5_df["SAMPLE"] == old_sample) & (fast5_df["RUN_ID"] == run_id))]
fast5_df.to_csv("ont_fast5_table.tsv", sep='\t', index=False)

bc_df = pd.read_csv("ont_basecall.tsv", sep='\t')
bc_df = bc_df.loc[~((bc_df["SAMPLE"] == old_sample) & (bc_df["RUN_ID"] == run_id))]
bc_df.to_csv("ont_basecall.tsv", sep='\t', index=False)

copy_rename = os.path.join(new_prefix, "/".join(args.input.split("/")[1:]))

os.makedirs(os.path.dirname(copy_rename), exist_ok=True)

df[out_cols].to_csv(copy_rename, sep="\t", index=False)

with open(output_file, "w") as outfile:
    outfile.write("\n".join(migrated_files)+"\n")

logging.info(f"Output written to {output_file}")
