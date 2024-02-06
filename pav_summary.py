#!/bin/env python

import argparse
import logging
import pandas as pd
import numpy as np
import os


def find_files(args):
    file_dict = {}
    for sample in args.samples:
        for directory in args.pav_dir:
            if sample not in file_dict:
                if os.path.isfile(
                    os.path.join(directory, f"results/{sample}/bed", f"snv_snv.bed.gz")
                ):
                    file_dict[sample] = directory
                else:
                    continue
            else:
                continue
    return file_dict


def calc_stats(file_dict):
    all_df = pd.DataFrame()
    for sample in file_dict:
        desired_stats = {}
        desired_stats["SAMPLE"] = [sample]
        for file in ["snv_snv", "sv_del", "sv_ins", "sv_inv", "indel_del", "indel_ins"]:
            pav_file = os.path.join(
                file_dict[sample], f"results/{sample}/bed", f"{file}.bed.gz"
            )
            pav_df = pd.read_csv(pav_file, sep="\t")
            if file == "snv_snv":
                desired_stats["SNV_COUNT"] = [len(pav_df)]
            else:
                desired_stats[f"{file.upper()}_COUNT"] = [len(pav_df)]
                desired_stats[f"{file.upper()}_MEAN"] = [np.mean(pav_df["SVLEN"])]
                desired_stats[f"{file.upper()}_MEDIAN"] = [np.median(pav_df["SVLEN"])]
                desired_stats[f"{file.upper()}_TBP"] = [np.sum(pav_df["SVLEN"])]
                desired_stats[f"{file.upper()}_MAX"] = [np.max(pav_df["SVLEN"])]
                desired_stats[f"{file.upper()}_MIN"] = [np.min(pav_df["SVLEN"])]
        all_df = pd.concat([all_df, pd.DataFrame.from_dict(desired_stats)])
    return all_df


def pretty(all_df):
    column_dict = {
        "indel_del": "Indel (DEL)",
        "indel_ins": "Indel (INS)",
        "sv_del": "SV (DEL)",
        "sv_ins": "SV (INS)",
    }
    pretty_df = pd.DataFrame()
    pretty_df["Sample"] = all_df["SAMPLE"]
    pretty_df["SNV Count"] = all_df["SNV_COUNT"].apply(
        lambda x: "{:.3f}M".format(x / 1e6)
    )
    for column in column_dict:
        label = column_dict[column]
        pretty_df[f"{label} Count"] = all_df[f"{column.upper()}_COUNT"].apply(
            lambda x: "{:,.0f}".format(x)
        )
        pretty_df[f"{label} Med|Mean (bp)"] = all_df.apply(
            lambda row: "|".join(
                [
                    "{:,.0f}".format(row[x])
                    for x in [f"{column.upper()}_MEDIAN", f"{column.upper()}_MEAN"]
                ]
            ),
            axis=1,
        )
        if "sv_" in column:
            pretty_df[f"{label} Total len (Mbp)"] = all_df[f"{column.upper()}_TBP"].apply(
                lambda x: "{:.3f}".format(x / 1e6)
            )
    return pretty_df


def main(args):
    # Your main logic goes here
    logging.debug("Started processing with arguments: %s", args)
    file_dict = find_files(args)
    all_df = calc_stats(file_dict)
    all_df.to_csv(args.output, sep="\t", index=False)
    if args.pretty_table:
        pretty_df = pretty(all_df)
        pretty_df.to_csv(args.pretty_table, sep="\t", index=False)


if __name__ == "__main__":
    # Set up argparse for handling command line arguments
    parser = argparse.ArgumentParser(description="Your script description")

    # Example argument with multiple values
    parser.add_argument(
        "-d",
        "--pav_dir",
        nargs="+",
        type=str,
        help="Path to find PAV calls (can be multiple)",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--samples",
        nargs="+",
        type=str,
        help="Sample names to look for",
        required=True,
    )
    parser.add_argument("-o", "--output", type=str, help="Output file", required=True)
    parser.add_argument(
        "-p", "--pretty_table", type=str, help="Output file", required=False
    )

    # Set logging configuration
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Parse command line arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args)
