#!/usr/bin/env python3
"""
Usage: ./get_pb_stats.py --sample CHM1 --cohort pop --prefix /path/to/LRA
Author: Mei Wu, https://github.com/projectoriented

This is a script tuned for CCS/HiFi or Revio generated .fastq.gz files.
"""

import os
import glob
import sys
import logging
import argparse

import pandas as pd
import numpy as np

LOG = logging.getLogger()
logging.basicConfig(stream=sys.stdout, level="INFO", format='%(asctime)s - %(levelname)s - %(message)s')


def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    required = parser.add_argument_group('required')

    required.add_argument('--sample', required=True, help='Sample')
    required.add_argument('--cohort', required=True, help='Cohort, e.g. clinical, pop')
    required.add_argument('--prefix', required=True, help='Common directory path of all cohorts in the LRA')

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    sample = args.sample
    cohort = args.cohort
    prefix = args.prefix
    pb_obj = FindPB(prefix=prefix, sample=sample, cohort=cohort)

    pb_obj.apply_regex()
    pb_obj.make_directory()

    outpath = pb_obj.outpath

    final_df = pd.DataFrame()

    for idx, row in pb_obj.df.iterrows():
        final_df = pd.concat(
            [
                final_df,
                CalculateStats(filepath=[row.filepath], cell_name=row.filename).stats
            ]
        )

    # Add a total row
    fastqs = pb_obj.df.filepath.tolist()
    total_df = CalculateStats(filepath=fastqs, cell_name="total").stats
    final_df = pd.concat([final_df, total_df])

    del pb_obj
    final_df.to_csv(outpath, sep='\t', index=False, header=True)
    LOG.info(f"Wrote: {outpath}")


class CalculateStats:
    def __init__(self, filepath: list, cell_name: str):
        self.filepath = filepath
        self.cell_name = cell_name
        self.genome = 3.1

    @staticmethod
    def get_n50(vals):
        """copy/pasta: https://github.com/EichlerLab/ccs_stats/blob/master/rules/ccs_stats.snakefile"""
        vals = vals.sort_values(ascending=False)
        vals_csum = np.cumsum(vals)
        return vals.iloc[np.sum(vals_csum <= (vals_csum.iloc[-1] // 2))] / 1000

    @property
    def df(self):
        working_df = pd.concat(
            [
                pd.read_csv(file + '.fai', sep='\t', header=None, usecols=[0, 1]) for file in self.filepath
            ]
        )
        return working_df

    @property
    def stats(self) -> pd.DataFrame:
        working_df = self.df

        len_list = pd.Series(working_df[1].copy())
        len_list.sort_values(ascending=False, inplace=True)

        len_list_k = {
            "15k": pd.Series(working_df.loc[working_df[1] >= 15000][1].copy()),
            "18k": pd.Series(working_df.loc[working_df[1] >= 18000][1].copy())
        }

        coverage = np.sum(len_list) / (self.genome * 1000000000)
        coverage_k = {
            "15k": np.sum(len_list_k["15k"]) / (self.genome * 1000000000),
            "18k": np.sum(len_list_k["18k"]) / (self.genome * 1000000000),
        }

        return pd.DataFrame(data={
            'CELL': [self.cell_name],
            'COVERAGE_X': ["{:.2f}".format(coverage)],
            'COVERAGE_X_15_Kbp': ["{:.2f}".format(coverage_k["15k"])],
            'COVERAGE_X_18_Kbp': ["{:.2f}".format(coverage_k["18k"])],
            'READS': [len(len_list)],
            'N50_Kbp': ["{:.2f}".format(self.get_n50(len_list))],
            'FILES': ','.join(self.filepath)
        })


class FindPB:
    regex = r'(?P<common_dir>.*PacBio_HiFi)/(?P<filename>.*fastq.gz)'

    def __init__(self, prefix, sample, cohort):
        self.prefix = prefix
        self.sample = sample
        self.cohort = cohort
        self.glob_list = glob.glob(
            f"{os.path.join(prefix, cohort, sample)}/raw_data/PacBio_HiFi/*.fastq.gz")
        self.df = pd.DataFrame(data=self.glob_list, columns=["filepath"])

    @property
    def unique_indicies(self):
        return self.df.index.unique().tolist()

    def apply_regex(self):
        # Check if the DF is empty.
        if self.df.empty:
            raise OSError(f"Cannot find any PacBio files for {os.path.join(self.prefix, self.cohort, self.sample)}")

        regex_df = self.df.filepath.str.extract(self.__class__.regex, expand=True)
        self.df = pd.concat([self.df, regex_df], axis=1)

        self.df["directory_make"] = self.df.apply(lambda row: os.path.join(row.common_dir, "quick_stats"), axis=1)
        self.df["n50_filename"] = self.df.apply(lambda row: f"n50_{self.sample}.tsv", axis=1)

        # re-do the filename
        self.df["filename"] = self.df["filename"].apply(lambda x: x.split(".")[0])

        self.df.set_index("n50_filename", inplace=True)

        return self.df

    def make_directory(self):
        list_of_dirs = self.df.directory_make.unique().tolist()
        for v in list_of_dirs:

            # Check if we can write in the directory
            one_level_up = os.path.dirname(v)
            if not os.access(one_level_up, os.W_OK):
                raise RuntimeError(f"Cannot write in {v}")

            if not os.path.exists(v):
                LOG.info(f"Making {v} directory")
                os.makedirs(v)
                os.chmod(v, 0o775)
            else:
                LOG.info(f"{v} exists, skipping.")

    @property
    def outpath(self):
        if len(self.unique_indicies) > 1:
            raise OSError(f"There are multiple output names: {self.unique_indicies}")

        n50_filename = self.unique_indicies.pop()
        n50_abs_path = os.path.join(self.df.directory_make[0], n50_filename)
        return n50_abs_path


if __name__ == "__main__":
    sys.exit(main())
