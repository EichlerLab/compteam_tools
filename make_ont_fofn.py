#!/usr/bin/env python3
"""
Usage: ./make_ont_fofn.py --sample 14455_p1 --proj_dir absolute/path/to/clinical --output 14455_p1.fofn --filter_string 'lib=STD;model=sup-prom;bc=guppy;ver=6;ftype=bam'
Author: Mei Wu, https://github.com/projectoriented

This script will get you all the fastqs that belong to guppy ver 6.x.x for STD library
Slack/email me (Mei) if you experience any problems or want additional features
"""

from glob import iglob
import re
import pandas as pd
import sys
import argparse
import os


def main():
    """Run script to output ont fofn"""
    parser = get_parser()
    args = parser.parse_args()

    make_ont_fofn(sample=args.sample, fn=args.output, prefix=args.proj_dir, fltr_str=args.filter_string)


def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    required = parser.add_argument_group('required')

    required.add_argument('--sample', type=str, required=True,
                          help='e.g. 14455_p1')
    required.add_argument('--proj_dir', type=str, required=True,
                          help='Absolute path for project directory.')
    parser.add_argument('--filter_string', type=str, required=False,
                        help='A filter string of k=v delimited by ; and have any of these keys: (bc, model, ver, lib)')
    parser.add_argument('--output', '-o', type=str, required=False, default=sys.stdout,
                        help='Output')

    return parser


def make_ont_fofn(sample, fn, prefix, fltr_str):
    fltr_dict = {
        'bc': 'guppy',
        'model': 'sup-prom',
        'ver': '6',
        'lib': 'STD',
        'ftype': 'fastq'
    }
    if fltr_str:
        alt_fltr_list = fltr_str.split(';')
        alt_fltr_dict = dict([x.split('=') for x in alt_fltr_list])
        fltr_dict.update(alt_fltr_dict)

    basecaller, model, version, library, file_type = tuple(fltr_dict.values())

    if file_type == 'fastq':
        file_type = 'fastq.gz'
    elif file_type == 'bam':
        file_type = 'bam'
    else:
        sys.exit(f'Invalid ftype {file_type}. Choose from (fastq, bam)')

    regex = fr'(.*nanopore)/(?P<lib>{library})/(.*fastq)/(?P<runid>.*)/(?P<bc>{basecaller})/(?P<ver>{version}.*)/(?P<model>{model})/(.*fastq_pass)(?P<ftype>.+)'
    search_pattern = re.compile(regex)
    a = []

    search_path = os.path.join(prefix, sample, f'raw_data/nanopore/*/**/*fastq_pass.{file_type}')

    file_list = iglob(search_path, recursive=True)
    for f in file_list:
        match = search_pattern.match(f)
        if match:
            match_dict = match.groupdict()
            match_dict['fpath'] = f
            a.append(match_dict)

    fofn_df = pd.DataFrame(a)

    if fofn_df.empty:
        sys.exit(f'No {file_type} matched your filters for {sample}- please try again :)')

    # Make sure we follow the real path and not soft links.
    fofn_df["is_link"] = fofn_df.apply(lambda row: 1 if os.path.islink(row.fpath) else 0, axis=1)
    fofn_df = fofn_df.query("is_link == 0").reset_index(drop=True).drop(columns=["is_link"])

    # Make sure we have the real path from hard links.
    for row in fofn_df.itertuples():
        try:
            if os.stat(row.fpath).st_nlink > 1:
                fofn_df.loc[row.Index, "true_path"] = os.path.realpath(row.fpath)

                fofn_df.loc[row.Index, "fpath"] = os.path.realpath(row.fpath)
            else:
                fofn_df.loc[row.Index, "true_path"] = row.fpath
        except FileNotFoundError:
            fofn_df.drop([row.Index], inplace=True)

    fofn_df = fofn_df.drop_duplicates(subset=["true_path"], keep="first").drop(columns=["true_path"])

    # Drop the duplicate samples but keeping the one run with latest version
    fofn_df.sort_values(['ver'], inplace=True)
    fofn_df.drop_duplicates(['lib', 'runid', 'bc', 'model'], inplace=True, keep='last')

    return fofn_df['fpath'].to_csv(fn, header=False, index=False)


if __name__ == "__main__":
    sys.exit(main())
