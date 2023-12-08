import sys
from pathlib import Path
from collections import defaultdict
import subprocess

def main():
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(0)
    list_of_samples = sys.argv[1]

    with open(list_of_samples) as f:
        samples_in = [z.strip() for z in f.readlines()]

    qstat_cmd = (
        'for j in $(qstat -u $USER | tail -n +3 | cut -f 2 -d " "); '
        'do qstat -j $j; done'
    )
    qstat_all_str = str(subprocess.run(
        qstat_cmd, shell=True, check=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    ).stdout)
    qstat_jobs = [
        parse_qstat_str(z) for z in qstat_all_str.strip().split(
        "=============================================================="
        )
    ]
    qstat_jobs = [z for z in qstat_jobs if z]

    cromwell_dirs_per_sample = defaultdict(set)
    for p in Path('.').glob(
        'cromwell-executions/FlaggerEndToEnd/*/call-project/'
        'runProjectBlocksForFlagger/*/call-bam2pafHap1/inputs/*/*'
    ):
        cromwell_dir = str(p.parts[2])
        sample = p.stem
        if sample in samples_in:
            cromwell_dirs_per_sample[sample].add(cromwell_dir)

    for sample in samples_in:
        print_sample_status(
            sample, cromwell_dirs_per_sample[sample], qstat_jobs
        )


def print_usage():
    print(
        '''------
Usage: python monitor_flagger_jobs.py <sample_file>

Run this from Flagger base directory (containing cromwell-executions)
to monitor Flagger jobs running from there. Requires python version >=3.4.

sample_file: File with one sample name per line to check, no header. E.g.
  HG01457-hifiasm
  HG02011-hifiasm
------'''
)


def print_sample_status(sample, cromwell_dirs, qstat_jobs):
    end_to_end_jobs = []
    subjobs = []
    for j in qstat_jobs:
        if j['job_name'] == f'run_flagger_end_to_end_{sample}':
            end_to_end_jobs.append(j)
        elif sample in j['submit_cmd']:
            subjobs.append(j)

    # Process cromwell-executions dirs
    final_output_files = []
    for cromwell_dir in cromwell_dirs:
        for final_output_file in (
            'cromwell-executions/FlaggerEndToEnd'/Path(cromwell_dir)
        ).glob(f'**/output/{sample}*.alt_removed.flagger_final.bed'):
            final_output_files.append(str(final_output_file))
        subjobs += [
            z for z in qstat_jobs if (
                (cromwell_dir in z['stdout_path_list']) or 
                (cromwell_dir in z['submit_cmd'])
            )
        ]

    print (f'{sample}:\n')
    print(
        'cromwell-executions directories: ' + (', '.join(cromwell_dirs) if (
            cromwell_dirs
        ) else 'None')
    )
    print(
        'Final output files: ' + (', '.join(final_output_files) if (
            final_output_files
        ) else 'None')
    )

    if end_to_end_jobs:
        print(
            'End-to-end job: ' + ', '.join([
                f"{z['job_name']}: Job {z['job_number']}- state {z['job_state']}"
                for z in end_to_end_jobs
            ])
        )
    if subjobs:
        print(
            'Subjobs:\n' + '\n'.join([
                f"{z['job_name']}: Job {z['job_number']}- state {z['job_state']}"
                for z in subjobs
            ])
        )
    if end_to_end_jobs and subjobs:
        print('Summary: Running')
    elif end_to_end_jobs:
        print('Summary: Warning: End-to-end job running but no subjobs running')
    elif subjobs:
        print('Summary: Warning: Subjobs running but no end-to-end jobs running')
    elif final_output_files:
        print('Summary: Done')
    else:
        print('Summary: Failed')
    print('\n------\n')


def parse_qstat_str(qstat_str):
    qstat_str = qstat_str.split('\\n')
    return {
        z[:z.index(':')].split()[0].strip() : z[z.index(':') + 1:].strip()
        for z in qstat_str if ':' in z
    }


if __name__ == '__main__':
    main()