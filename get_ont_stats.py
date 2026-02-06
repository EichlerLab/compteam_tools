#!/usr/bin/env python3
"""
Usage: ./get_ont_stats.py --sample CHM1 --cohort pop --prefix /path/to/LRA
Author: Mei Wu, https://github.com/projectoriented
Modified: Youngjun Kwon
"""

import os
import glob
import sys
import logging
import argparse
from collections import defaultdict

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
    ont_obj = FindONT(prefix=prefix, sample=sample, cohort=cohort)

    ont_obj.apply_regex()
    ont_obj.make_directory()

    fastq_dict = ont_obj.get_fastqs()

    for k in fastq_dict.keys():
        final_df = pd.DataFrame()
        working_df = fastq_dict[k]
        if isinstance(working_df, pd.Series):
            working_df = working_df.to_frame().T

        for idx, row in working_df.iterrows():
            final_df = pd.concat(
                [
                    final_df,
                    CalculateStats(filepath=[row.filepath], cell_name=row.run_id).stats
                ]
            )
        # Add a total row
        fastqs = working_df.filepath.tolist()
        total_df = CalculateStats(filepath=fastqs, cell_name="total").stats
        final_df = pd.concat([final_df, total_df])

        final_df.to_csv(k, sep='\t', index=False, header=True)
        LOG.info(f"Wrote: {k}")


class CalculateStats:
    def __init__(self, filepath: list, cell_name: str):
        self.filepath = filepath
        self.cell_name = cell_name
        self.genome = 3.1

    @staticmethod
    def get_n50(vals):
        """copy/pasta https://github.com/EichlerLab/compteam_tools/blob/main/ont_stats"""
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
        """copy/pasta https://github.com/EichlerLab/compteam_tools/blob/main/ont_stats"""
        working_df = self.df

        len_list = pd.Series(working_df[1].copy())
        len_list.sort_values(ascending=False, inplace=True)

        len_list_k = pd.Series(working_df.loc[working_df[1] >= 100000][1].copy())

        coverage = np.sum(len_list) / (self.genome * 1000000000)
        coverage_k = np.sum(len_list_k) / (self.genome * 1000000000)

        return pd.DataFrame(data={
            'CELL': [self.cell_name],
            'COVERAGE_X': ["{:.2f}".format(coverage)],
            'COVERAGE_X_100_Kbp': ["{:.2f}".format(coverage_k)],
            'READS': [len(len_list)],
            'N50_Kbp': ["{:.2f}".format(self.get_n50(len_list))],
            'FILES': ','.join(self.filepath)
        })


class FindONT:
    regex = r'(?P<common_dir>.*nanopore)/(?P<library>[A-Z]{2,5})/fastq/(?P<run_id>.+)/(?P<basecaller>.+)/(?P<version>\d+.\d+.\d+)/(?P<model>.+)/(?P<filename>.*_pass.*fastq.gz)'

    def __init__(self, prefix, sample, cohort):
        self.prefix = prefix
        self.sample = sample
        self.cohort = cohort
        # added a condition that fastq have filtered and has .filt. in the name.
        glob_list_tmp = glob.glob(
            f"{os.path.join(prefix, cohort, sample)}/raw_data/nanopore/*/fastq/*/*/*/*/*_pass*fastq.gz")

        groups = defaultdict(list)

        for f in glob_list_tmp:
            if f.endswith("pass.filt.fastq.gz"):
                key = f.replace("pass.filt.fastq.gz", "")
            elif f.endswith("pass.fastq.gz"):
                key = f.replace("pass.fastq.gz", "")
            else:
                continue
            groups[key].append(f)

        glob_list = []
        for files in groups.values():
            filt = [file for file in files if file.endswith("pass.filt.fastq.gz")]
            if filt:
                glob_list.extend(filt)
            else:
                glob_list.extend(files)

        self.glob_list = [ fastq for fastq in glob_list if not "HERRO" in fastq ]
        self.df = pd.DataFrame(data=self.glob_list, columns=["filepath"])

    @property
    def unique_indicies(self):
        return self.df.index.unique().tolist()

    def apply_regex(self):
        # Check if the DF is empty.
        if self.df.empty:
            raise OSError(f"Cannot find any files for {os.path.join(self.prefix, self.cohort, self.sample)}")

        regex_df = self.df.filepath.str.extract(self.__class__.regex, expand=True)
        self.df = pd.concat([self.df, regex_df], axis=1)

        # Just take the major version number
        self.df["version"] = self.df["version"].str.split(".", n=1, expand=True)[0]

        self.df["directory_make"] = self.df.apply(
            lambda row: os.path.join(row.common_dir, row.library, "quick_stats"),
            axis=1)
        self.df["n50_filename"] = self.df.apply(lambda row: f"n50_{row.basecaller}_{row.model}_v{row.version}.txt",
                                                axis=1)

        self.df.set_index(["version", "model", "library", "directory_make", "n50_filename"], inplace=True, drop=False)
        self.df.sort_index(inplace=True)
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

    def get_fastqs(self):
        fastq_dict = {}
        for v in self.unique_indicies:
            dict_key = os.path.join(v[3], v[-1])  # directory_make + n50_filename
            fastq_dict[dict_key] = self.df.loc[v, ["filepath", "run_id"]]
        return fastq_dict


if __name__ == "__main__":
    sys.exit(main())
