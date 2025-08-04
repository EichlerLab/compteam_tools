#!/bin/env python

import pandas as pd
import numpy as np
import sys
import os
import glob
import re
import argparse

LRA = "/net/eichler/vol28/projects/long_read_archive/nobackups"

def check_copy_record(cohort, sample, seqtype, runid):
    copy_record = f"{LRA}/{cohort}/{sample}/raw_data/nanopore/{seqtype}/copy_record/fast5/{runid}.tab.gz"
    if os.path.isfile(copy_record):
        return True
    else:
        return np.nan
    
def expand_to_rows(df, column_name):
    rows = []
    for _, row in df.iterrows():
        run_ids = row[column_name]
        for run_id in run_ids:
            new_row = row.copy()
            new_row[column_name] = run_id
            rows.append(new_row)
    return pd.DataFrame(rows)
    
def get_basecalled_bams(cohort, sample, seqtype, runid):
    glob_str = f"{LRA}/{cohort}/{sample}/raw_data/nanopore/{seqtype}/fastq/{runid}/*/*/*/*pass*bam"
    basecalled_bams = list(glob.iglob(glob_str))
    if len(basecalled_bams) == 0:
        return [np.nan]
    else:
        return basecalled_bams


def check_raw_data(cohort, sample, seqtype, runid):
    glob_str = f"{LRA}/{cohort}/{sample}/raw_data/nanopore/{seqtype}/fast5/{runid}/*/*5"
    try:
        raw_data = list(glob.iglob(glob_str))[0]
        return os.path.basename(raw_data).split(".")[-1].upper()
    except:
        return "REMOVED"



def get_basecall_info(basecalled_bam):
    try:
        filename = os.path.basename(basecalled_bam)
        basecaller, version, model_dir = basecalled_bam.split("/")[14:17]        
    except:
        return "NA", "NA", "NA"
    if basecaller == "guppy":
        if "kit14" in model_dir: # R10
            if "sup" in model_dir:
                model = "dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_sup_prom"
            elif "hac" in model_dir:
                model = "dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac_prom"
        else: # R9
            if "sup" in model_dir:
                model = "dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_sup_prom"
            elif "hac" in model_dir:
                model = "dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_hac_prom"
    elif basecaller == "dorado":
        if float(version.replace(".","")) <= float("0.4.2".replace(".","")): # dorado 0.4.2
            if "kit14" in "sup-prom": # R10
                model = "dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG"
            else: # R9
                model = "dna_r9.4.1_e8_sup@v3.3_5mCG_5hmCG"
        else: # dorado version > 0.4.2
            model = filename.split(f"{basecaller}-{version}-")[1].replace("_pass.bam","")
    return basecaller, version, model
    

def check_fastq(basecalled_bam):

    try:
        expected_fastq = basecalled_bam.replace(".bam", ".fastq.gz")
    except:
        return False, False

    if "_fastq_pass." in expected_fastq:
        alt_fastq = expected_fastq.replace("_fastq_pass.", "_pass.")
    else:
        alt_fastq = expected_fastq.replace("_pass.", "_fastq_pass.")

    for fastq in [expected_fastq, alt_fastq]:
        if os.path.isfile(fastq):
            exist = True
            complete = True if os.path.isfile(f"{fastq}.fai") else False
            return exist, complete
        
    return False, False


def get_chemisty(profile):
    if re.search("kit14", profile):
        return "R10"
    else:
        return "R9"

def main():
    if sys.stdin.isatty():
        if args.infile is None:
            parser.error("Input file not provided and no pipe input detected.")
        df = pd.read_csv(args.infile, sep="\t", header=None, names = ["SAMPLE","SEQTYPE","RUNID","PROFILE"])
    else:
        df = pd.read_csv(sys.stdin, sep="\t", header=None, names = ["SAMPLE","SEQTYPE","RUNID","PROFILE"])
    
    df["COHORT"] = args.cohort
    df["CHEMISTRY"] = df["PROFILE"].apply(get_chemisty)
    df["COPY_RECORD_EXISTS"] = df.apply(lambda row:check_copy_record(row["COHORT"], row["SAMPLE"],  row["SEQTYPE"], row["RUNID"]),axis=1)
    df["BASECALLED_BAM"] = df.apply(lambda row:get_basecalled_bams(row["COHORT"], row["SAMPLE"],  row["SEQTYPE"], row["RUNID"]),axis=1)
    df["RAW_DATA"] = df.apply(lambda row:check_raw_data(row["COHORT"], row["SAMPLE"],  row["SEQTYPE"], row["RUNID"]),axis=1)
    df = expand_to_rows(df.copy(), "BASECALLED_BAM")
    df[["BASECALLER", "VERSION", "MODEL"]] = df["BASECALLED_BAM"].apply(get_basecall_info).apply(pd.Series)
    df[["FASTQ_EXIST", "FASTQ_COMPLETE"]] = df["BASECALLED_BAM"].apply(lambda x: pd.Series(check_fastq(x), dtype=bool))

    df = df[["COHORT","SAMPLE","CHEMISTRY","SEQTYPE","RUNID","BASECALLER","VERSION","MODEL","FASTQ_EXIST","FASTQ_COMPLETE","RAW_DATA"]]

    df.to_csv(sys.stdout, sep="\t", index=False)
    


    # df["BASECALL"] = df.apply(lambda row:check_basecall(row["COHORT"], row["LRA_SAMPLE"], row["SEQTYPE"], row["FC_ID"]),axis=1)
    # print (df)
    # df.to_csv(outname, sep="\t", index=False)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True, help='Cohort name [clinical | pop | nhp | external]', choices=['clinical', 'pop', 'nhp', 'external'])
    parser.add_argument('infile', nargs='?', help='Input TSV (ont_basecall.tsv) file (optional, or use pipe input)')
    args = parser.parse_args()
    main()




        
    # check_dict = dict()
    # for fastq in basecalled_fastqs:
    #     fastq_basename = os.path.basename(fastq)
    #     basecaller = fastq.split("/")[14]
    #     version = fastq.split("/")[15]
    #     model_dir = fastq.split("/")[16]
    #     check_dict.setdefault(basecaller, dict())

    #     try:
    #         model_str = re.findall('_[a-z_]+_v[0-9.]+',fastq_basename)[0]
    #         model = model_str.split("_v")[0][1:]
    #         model_verison = re.findall('v[0-9.]+',model_str)[0]
    #         modbase = "_".join((fastq_basename.split(model_str)[1]).split("-")[1].split("_")[:-1])
    #         full_model_info = f"{model}_{model_verison}-{modbase}"
    #     except:
    #         if basecaller == "guppy":
    #             if "kit14" in model_dir: # R10
    #                 if "sup" in model_dir:
    #                     full_model_info = "r10_sup-5mCG_5hmCG"
    #                 elif "hac" in model_dir:
    #                     full_model_info = "r10_hac-5mCG_5hmCG"
    #             else: # R9
    #                 if "sup" in model_dir:
    #                     full_model_info = "r9_sup-5mCG_5hmCG"
    #                 elif "hac" in model_dir:
    #                     full_model_info = "r9_hac-5mCG_5hmCG"
    #         else: # dorado <= 0.4.2
    #             if model_dir == "sup-prom": # R9
    #                 full_model_info = "r9_sup_v3.3-5mCG_5hmCG"
    #             else:
    #                 full_model_info = "r10_sup_v4.2.0-5mCG_5hmCG"
        
    #     check_dict[basecaller].setdefault(version, [])
    #     check_dict[basecaller][version].append(full_model_info)
    
    # basecall_info_all = []
    # for basecaller in sorted(check_dict.keys()):
    #     for version in sorted(check_dict[basecaller].keys()):
    #         model_info = ",".join(list(set(check_dict[basecaller][version])))
    #         basecall_info_all.append(f"{basecaller}-{version}:{model_info}")
    # basecall_info = "|".join(basecall_info_all)
    # return basecall_info    



