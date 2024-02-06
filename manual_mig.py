#!/bin/env python

from datetime import datetime
import collections
import hashlib
import pandas as pd
import numpy as np
import re
import shutil
import stat
import gzip
import io
import glob
import xml.etree.ElementTree as ET
import math
import logging
import os, sys
import argparse


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


parser = argparse.ArgumentParser(
    description="Python script to manually migrate CCS data from an Excel file",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)


parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="Input EXCEL file",
)

args = parser.parse_args()


logging.basicConfig(
    format="%(levelname)s (%(asctime)s): %(message)s (Line: %(lineno)d [%(filename)s])",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.WARNING,
)

# args = pd.Series(["231009_manualLRAmoves.xlsx"], index=["input"])

if args.input.endswith("xlsx"):
    df = pd.read_excel(args.input)
elif args.input.endswith("tab") or args.input.endswith("tsv"):
    df = pd.read_csv(args.input, sep="\t")
else:
    logging.error("Unrecognized file type: file must end in xlsx, tsv, or tab")
    sys.exit(1)

try:
    df[["LRA_Sample", "Path", "LRA_dir"]]
except KeyError:
    logging.error(
        "LRA_Sample, Path, LRA_dir are required as columns for the input data"
    )
    sys.exit(1)

df["Type"] = df["Path"].apply(lambda val: val.split(".")[-1])
df = df.loc[df["Type"].str.contains("bam")]
df["BC"] = df["Path"].apply(lambda val: val.split(".")[-2].split("-")[0])

df_xml = df.loc[df["Type"] == "bam"].copy()

df_xml["Type"] = "xml"
df_xml["Path"] = df_xml["Path"].str.replace(
    ".bam", ".consensusreadset.xml", regex=False
)

df_pbi = df.loc[df["Type"] == "bam"].copy()

df_pbi["Type"] = "pbi"
df_pbi["Path"] = df_pbi["Path"] + ".pbi"


df_all = pd.concat([df, df_xml, df_pbi]).reset_index(drop=True)

df_all["file_stat"] = df_all["Path"].apply(
    lambda val: os.stat(val) if os.path.isfile(val) else "SKIP"
)


df_all = df_all.loc[df_all["file_stat"] != "SKIP"].copy()

df_all["RUN_ID"] = df_all["Path"].apply(
    lambda val: "_".join(os.path.basename(val).split("_")[0:3]).replace("m", "r")
)
df_all["CELL"] = "M01"
df_all["SEQ_TYPE"] = df_all["RUN_ID"].apply(
    lambda val: "revio" if val.startswith("r8") else "ccs"
)
df_all["CCS_JOB_ID"] = "MANUAL_DEMUX"
df_all["BARCODE"] = df_all["BC"]
df_all["SOURCE_PATH"] = df_all["Path"]
df_all["MD5"] = df_all["Path"].apply(lambda val: _get_checksum(val))
df_all["STATUS"] = "Copied"
df_all["SAMPLE"] = df_all["LRA_Sample"]
df_all["DEST_PATH"] = df_all.apply(
    lambda row: os.path.join(
        "/net/eichler/vol28/projects/long_read_archive/nobackups/",
        row["LRA_dir"],
        row["SAMPLE"],
        "raw_data/PacBio_HiFi",
        os.path.basename(row["Path"]),
    )
    if row["Type"] != "xml"
    else os.path.join(
        "/net/eichler/vol28/projects/long_read_archive/nobackups/",
        row["LRA_dir"],
        row["SAMPLE"],
        "raw_data/PacBio_HiFi",
        "reports",
        os.path.basename(row["Path"]),
    ),
    axis=1,
)
df_all["SIZE"] = df_all["file_stat"].apply(lambda val: val.st_size)
df_all["MOD_TIME"] = df_all["file_stat"].apply(
    lambda val: datetime.fromtimestamp(val.st_mtime).strftime("%Y-%m-%d %H:%M:%S.%f")
)


columns_out = [
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
    "MD5",
]

df_all = df_all.sort_values(["LRA_Sample", "BC"])
df_all = df_all.set_index(["LRA_Sample", "BC"])

for index in df_all.index:
    df_copy_record = df_all.loc[index].copy().reset_index()
    for index_copy in df_copy_record.index:
        os.makedirs(
            os.path.dirname(df_copy_record.at[index_copy, "DEST_PATH"]), exist_ok=True
        )
        print(df_copy_record.at[index_copy, "DEST_PATH"])
        shutil.copy2(
            df_copy_record.at[index_copy, "SOURCE_PATH"],
            df_copy_record.at[index_copy, "DEST_PATH"],
        )
    sample = df_copy_record.iloc[0]["SAMPLE"]
    barcode = "-" + df_copy_record.iloc[0]["BARCODE"]
    cohort = df_copy_record.iloc[0]["LRA_dir"]
    seq_type = df_copy_record.iloc[0]["SEQ_TYPE"]
    run_id = df_copy_record.iloc[0]["RUN_ID"]
    cell_id = df_copy_record.iloc[0]["CELL"]
    RECORD_TAB_PATTERN = f"/net/eichler/vol28/projects/long_read_archive/nobackups/{cohort}/{sample}/raw_data/PacBio_HiFi/copy_record/{seq_type}/{run_id}/{cell_id}{barcode}.tab.gz"
    if os.path.isfile(RECORD_TAB_PATTERN):
        df_out = df_copy_record[columns_out].append(pd.read_csv(RECORD_TAB_PATTERN, sep='\t'))
    else:
        df_out = df_copy_record.copy()
    os.makedirs(os.path.dirname(RECORD_TAB_PATTERN), exist_ok=True)
    df_out[columns_out].to_csv(RECORD_TAB_PATTERN, sep="\t", index=False)
    print(RECORD_TAB_PATTERN)

