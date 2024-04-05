import sys
import os
import shutil
import argparse
from pathlib import Path
import csv


def main():
    args = parse_args(sys.argv[1:])

    validate_manifest(args.manifest, args.no_header, args.add_bais)

    copy_files_to_LRA(
        args.manifest, args.no_header, args.lra_target_path, args.add_bais,
        args.dry_run
    )


def parse_args(args):
    p = argparse.ArgumentParser(
        description='Copy filepaths defined in a manifest to specified long '
        'read archive path.'
    )
    p.add_argument(
        'manifest', type=Path,
        help='Path to a manifest TSV with first column sample and second '
        'column filepath to migrate. Assumes header present with any values.'
    )
    p.add_argument(
        'lra_target_path', type=str,
        help='Path to copy the files to, where "{sample}" will be replaced '
        'with the sample name. Example: /net/eichler/vol28/projects/'
        'long_read_archive/nobackups/clinical/{sample}/alignments/Illumina/'
        'CHM13_v2.0/alignment.bam'
    )
    p.add_argument(
        '--add_bais', '-a', action='store_true',
        help='Also copy over the .bai index file for all files in the '
        'manifest, assuming bams.'
    )
    p.add_argument(
        '--dry_run', '-d', action='store_true',
        help='Print each source and target filepath to be copied, without '
        'actually copying over the files.'
    )
    p.add_argument(
        '--no_header', '-n', action='store_true', help='Manifest has no header.'
    )
    return p.parse_args()


def validate_manifest(manifest, no_header, add_bais):
    # Make sure paths to migrate exist, each sample only occurs
    # once, and filename includes sample name
    samples = set()
    with manifest.open() as f:
        r = csv.DictReader(f, delimiter='\t', fieldnames=['sample', 'fpath'])
        if not no_header:
            next(r)
        for line in r:
            fpath = Path(line['fpath'])
            sample = line['sample']
            if not fpath.exists():
                raise ValueError(f'Missing filepath in manifest: {fpath}')
            if add_bais:
                if not Path(str(fpath) + '.bai').exists():
                    raise ValueError(
                        f'Missing filepath in manifest: {fpath}.bai'
                    )
            if sample in samples:
                raise ValueError(f'Sample {sample} repeated in manifest')
            if sample not in str(fpath):
                raise ValueError(f'Sample {sample} not in filepath {fpath}')
            samples.add(sample)


def copy_files_to_LRA(manifest, no_header, lra_target_path, add_bais, dry_run):
    if dry_run:
        print('Dry run:')
    with manifest.open() as f:
        r = csv.DictReader(f, delimiter='\t', fieldnames=['sample', 'fpath'])
        if not no_header:
            next(r)
        for line in r:
            fpath = Path(line['fpath'])
            sample = line['sample']
            out_fpath = Path(lra_target_path.replace('{sample}', sample))
            out_fpath_bai = Path(f'{out_fpath}.bai')
            if out_fpath.exists():
                raise ValueError(f'Target filepath {out_fpath} already exists')
            if add_bais and out_fpath_bai.exists():
                raise ValueError(f'Target filepath {out_fpath_bai} already exists')
            os.makedirs(out_fpath.parent, exist_ok=True)
            print(f'Copying {fpath} to {out_fpath}')
            if not dry_run:
                shutil.copyfile(fpath, out_fpath)
            if add_bais:
                print(f'Copying {fpath}.bai to {out_fpath_bai}')
                if not dry_run:
                    shutil.copyfile(f'{fpath}.bai', out_fpath_bai)


if __name__ == '__main__':
    main()
