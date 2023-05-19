#!/usr/bin/env python3
"""
Usage: ./get_ont_stats.py --sample CHM1 --cohort pop --prefix /path/to/LRA
Author: Mei Wu, https://github.com/projectoriented
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
    """Get parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True, help='Sample')
    parser.add_argument('--prefix', required=False, help='Common directory path of all cohorts in the LRA')
    parser.add_argument('--cohort', required=True, help='Cohort, e.g. clinical, pop')

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    sample = args.sample
    cohort = args.cohort
    prefix = args.prefix
    ont_obj = OntStats(prefix=prefix, sample=sample, cohort=cohort)

    ont_obj.apply_regex()

    ont_obj.make_directory()
    fastq_dict = ont_obj.get_fastqs()

    for k in fastq_dict.keys():
        CalculateStats(files=fastq_dict[k], sample=sample).write_out(outpath=k)


class CalculateStats:
    """copy/pasta https://github.com/EichlerLab/compteam_tools/blob/main/ont_stats"""
    def __init__(self, files: list, sample):
        self.files = files
        self.sample = sample
        self.genome = 3.1

    @staticmethod
    def get_n50(vals):
        vals = vals.sort_values(ascending=False)
        vals_csum = np.cumsum(vals)
        return vals.iloc[np.sum(vals_csum <= (vals_csum.iloc[-1] // 2))] / 1000

    @property
    def df(self):
        working_df = pd.concat([pd.read_csv(file + '.fai', sep='\t', header=None, usecols=[0, 1]) for file in self.files])
        return working_df

    @property
    def stats(self) -> dict:
        working_df = self.df

        len_list = pd.Series(working_df[1].copy())
        len_list.sort_values(ascending=False, inplace=True)

        len_list_k = pd.Series(working_df.loc[working_df[1] >= 100000][1].copy())

        coverage = np.sum(len_list) / (self.genome * 1000000000)
        coverage_k = np.sum(len_list_k) / (self.genome * 1000000000)

        return {
                'SAMPLE': [self.sample],
                'COVERAGE_X': ["{:.2f}".format(coverage)],
                'COVERAGE_X_100_Kbp': ["{:.2f}".format(coverage_k)],
                'READS': [len(len_list)],
                'N50_Kbp': ["{:.2f}".format(self.get_n50(len_list))],
                'FILES': ','.join(self.files)
            }

    def write_out(self, outpath):
        stats_dict = self.stats
        out_df = pd.DataFrame.from_dict(stats_dict)
        out_df.to_csv(outpath, sep='\t', index=False)
        LOG.info(f"Wrote: {outpath}")


class OntStats:
    regex = r'(?P<common_dir>.*nanopore)/(?P<library>[A-Z]{2,5})/fastq/.+/(?P<basecaller>guppy)/(?P<version>\d+.\d+.\d+)/(?P<model>.+)/(?P<filename>.*_fastq_pass.fastq.gz)'

    def __init__(self, prefix, sample, cohort):
        self.prefix = prefix
        self.sample = sample
        self.cohort = cohort
        self.glob_list = glob.glob(
            f"{os.path.join(prefix, cohort, sample)}/raw_data/nanopore/*/fastq/*/*/*/*/*_fastq_pass.fastq.gz")
        self.df = pd.DataFrame(data=self.glob_list, columns=["filepath"])

    @property
    def unique_indicies(self):
        return self.df.index.unique().tolist()

    def apply_regex(self):
        regex_df = self.df.filepath.str.extract(self.__class__.regex, expand=True)
        self.df = pd.concat([self.df, regex_df], axis=1)

        try:
            # Just take the major version number
            self.df["version"] = self.df["version"].str.split(".", n=1, expand=True)[0]
        except KeyError:
            LOG.warning(f"Something went wrong with the path: {self.df.filepath}")

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
                os.makedirs(v, mode=0o775)
            else:
                LOG.info(f"{v} exists, skipping.")

    def get_fastqs(self):
        fastq_dict = {}
        for v in self.unique_indicies:
            dict_key = os.path.join(v[3], v[-1]) # directory_make + n50_filename
            try:
                list_of_fastqs = self.df.loc[v, "filepath"].tolist()
            except AttributeError:
                list_of_fastqs = [self.df.loc[v, "filepath"]]
            fastq_dict[dict_key] = list_of_fastqs
        return fastq_dict


if __name__ == "__main__":
    sys.exit(main())
