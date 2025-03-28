# Make gene output html to be converted to image
import sys
import re
import math
import json
from collections import defaultdict, namedtuple
import argparse
import numpy as np
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from yattag import Doc

OutputLine = namedtuple('OutputLine', [
    'seq', 'num_repeats', 'num_haps', 'interrs_str', 'padding_i', 'seq_no_tags',
    'num_trgt_haps'
])
TrimmedSequence = namedtuple('TrimmedSequence', [
    'trimmed_seq', 'repeat_seq', 'seq_pre_repeat', 'seq_post_repeat'
])

REPEAT_UNIT_SHAPE = '●'
MIN_HAPS_SHOW_BARS = 5
MIN_SEQ_LEN = 5

def main():
    args = parse_args()
    assert args.out.is_dir()
    global config
    with open(args.gene_json) as f:
        config = json.load(f)

    male_samples = []
    if config['sex_per_sample_fname'] is not None:
        sex_per_sample = {}
        with open(config['sex_per_sample_fname']) as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                assert line[0] not in sex_per_sample
                sex_per_sample[line[0]] = line[1]
        male_samples = set([z for z in sex_per_sample if sex_per_sample[z] == 'Male'])

    print('--Reading in ASM sequences--')
    asm_seq_per_hap = read_in_asm(config['asm_sequences_fname'])
    asm_haps_per_repeat_seq, asm_rep_seq_tuple_per_repeat_seq = get_haps_per_repeat_seq(asm_seq_per_hap, male_samples)

    if config['trgt_sequences_fname'] is not None:
        print('--Reading in TRGT sequences--')
        trgt_seq_per_hap  = read_in_trgt(config['trgt_sequences_fname'])  # TRGT sequence in repeat region
        trgt_haps_per_repeat_seq, trgt_rep_seq_tuple_per_repeat_seq = get_haps_per_repeat_seq(trgt_seq_per_hap, male_samples)

    print('--Parsing inputs...--')
    print('--All sequences--')
    (sorted_out_lines, subs_used, asm_repeat_seqs) = process_sequences(
        asm_haps_per_repeat_seq.keys(), asm_haps_per_repeat_seq, asm_rep_seq_tuple_per_repeat_seq,
        sort_by=args.sort_by, min_haps=args.min_haps, min_repeats=args.min_repeats
    )
    if config['trgt_sequences_fname'] is not None:
        print('--TRGT only--')
        (sorted_out_lines_trgt, subs_used_trgt, trgt_repeat_seqs) = process_sequences(
            trgt_haps_per_repeat_seq.keys(), trgt_haps_per_repeat_seq, trgt_rep_seq_tuple_per_repeat_seq,
            seqs_to_exclude=asm_repeat_seqs, sort_by=args.sort_by, min_haps=args.min_haps, min_repeats=args.min_repeats
        )
        print('--ASM only--')
        (sorted_out_lines_asm, subs_used_asm, _) = process_sequences(
            asm_haps_per_repeat_seq.keys(), asm_haps_per_repeat_seq, asm_rep_seq_tuple_per_repeat_seq,
            seqs_to_exclude=trgt_repeat_seqs, sort_by=args.sort_by, min_haps=args.min_haps, min_repeats=args.min_repeats
        )
        #print(sorted(list(trgt_repeat_seqs), key=lambda u:len(u))[:10])
        print('--ASM + TRGT match--')
        (sorted_out_lines_match, subs_used_match, _) = process_sequences(
            asm_haps_per_repeat_seq.keys(), asm_haps_per_repeat_seq, asm_rep_seq_tuple_per_repeat_seq,
            seqs_to_exclude=asm_repeat_seqs.symmetric_difference(trgt_repeat_seqs),  # in one set but not the other
            trgt_haps=trgt_haps_per_repeat_seq, sort_by=args.sort_by, min_haps=args.min_haps, min_repeats=args.min_repeats
        )

    print('--Writing HTMLs--')
    if asm_repeat_seqs:
        write_html(args.out / f'{config['gene_name']}_viz.html', sorted_out_lines, subs_used)
    if config['trgt_sequences_fname'] is not None:
        write_html(args.out / f'{config['gene_name']}_viz_TRGT_only.html', sorted_out_lines_trgt, subs_used_trgt)#, divide_into=3)
        if asm_repeat_seqs:
            write_html(args.out / f'{config['gene_name']}_viz_asm_only.html', sorted_out_lines_asm, subs_used_asm)#, divide_into=2)
            write_html(args.out / f'{config['gene_name']}_viz_matches_TRGT.html', sorted_out_lines_match, subs_used_match)#, divide_into=3)

    print('--Statistics--')
    if config['trgt_sequences_fname'] is not None:
        n_asm_seqs_confirmed = 0
        n_asm_seqs_not_confirmed = 0
        n_single_hap_asm_seqs_confirmed = 0
        n_single_hap_asm_seqs_not_confirmed = 0
        for repeat_seq in asm_haps_per_repeat_seq:
            n_haps = len(set(asm_haps_per_repeat_seq[repeat_seq]))
            samples = set([z.split('_')[0] for z in asm_haps_per_repeat_seq[repeat_seq]])
            trgt_samples = set([z.split('_')[0] for z in trgt_haps_per_repeat_seq[repeat_seq]])
            if samples.intersection(trgt_samples):
                n_asm_seqs_confirmed += 1
                if n_haps == 1:
                    n_single_hap_asm_seqs_confirmed += 1
            else:
                n_asm_seqs_not_confirmed += 1
                if n_haps == 1:
                    n_single_hap_asm_seqs_not_confirmed += 1
        print(f'{n_asm_seqs_confirmed} assembly sequences confirmed by TRGT')
        print(f'{n_asm_seqs_not_confirmed} assembly sequences NOT confirmed by TRGT')
        print(f'{n_single_hap_asm_seqs_confirmed} single-haplotype assembly sequences confirmed by TRGT')
        print(f'{n_single_hap_asm_seqs_not_confirmed} single-haplotype assembly sequences NOT confirmed by TRGT')


def parse_args():
    p = argparse.ArgumentParser(
        description='Make a figure with repeat sequences encoding according to '
        'specified string substitutions.'
    )
    p.add_argument(
        'gene_json', type=Path, help='Path to JSON defining input filepaths and '
        'substitutions for gene.'
    )
    p.add_argument(
        '-s', '--sort_by', type=str, help='Sort by either descending repeat '
        'length or number of haplotypes with repeat sequence. (hap)',
        default='hap', choices=['repeat', 'hap']
    )
    p.add_argument(
        '-m', '--min_haps', type=int, help='Only retain repeat sequences found '
        'in at least this number of haplotypes. (0)', default=0
    )
    p.add_argument(
        '-r', '--min_repeats', type=int, help='Only retain sequences with at '
        'least this many repeat units. (0)', default=0
    )
    p.add_argument(
        '-o', '--out', type=Path, help='Output directory for plot HTMLs. (.)',
        default=Path('.')
    )
    if len(sys.argv) < 2:
        p.print_help()
        sys.exit(0)
    return p.parse_args()


def get_haps_per_repeat_seq(seq_per_hap, male_samples):
    haps_per_repeat_seq = defaultdict(set)  # repeat seq: samples (phased) it occurs in
    males_two_haps = set()
    males_two_different_haps = set()
    rep_seq_tuple_per_repeat_seq = {} # full sequence representative for repeat seq

    # Preprocessing
    for hap in seq_per_hap:
        seq_tuple = seq_per_hap[hap]
        assert (hap.endswith('_hap1') or hap.endswith('_hap2'))
        if (
            (hap.endswith('_hap1')) and
            (hap.replace('_hap1', '') in male_samples) and
            (hap.replace('_hap1', '_hap2') in seq_per_hap)
        ):
            sample = hap.replace('_hap1', '')
            if seq_tuple.repeat_seq != seq_per_hap[hap.replace('_hap1', '_hap2')].repeat_seq:  # repeat region
                #print(sample, '\n', seq_per_hap[hap], '\n', seq_per_hap[hap.replace('_hap1', '_hap2')])
                males_two_different_haps.add(sample)
            males_two_haps.add(sample)

    print(f'Skipping {len(males_two_different_haps)} males with two different haplotypes')
    print(f'Condensing {len(males_two_haps) - len(males_two_different_haps)} males with duplicate haplotypes')

    for hap in seq_per_hap:
        sample = hap.split('_')[0]
        # males w >1 haplotype on chr X
        if not (
            ((hap.endswith('_hap2')) and (sample in males_two_haps)) or
            (sample in males_two_different_haps)
        ):
            seq_tuple = seq_per_hap[hap] # parsed sequence tuple for haplotype
            haps_per_repeat_seq[seq_tuple.repeat_seq].add(hap)
            #if seq_tuple.repeat_seq in rep_seq_tuple_per_repeat_seq and rep_seq_tuple_per_repeat_seq[seq_tuple.repeat_seq] != seq_tuple:
            #    print(f'Multiple sequences for {seq_tuple.repeat_seq}: {rep_seq_tuple_per_repeat_seq[seq_tuple.repeat_seq]}, {seq_tuple}')
            rep_seq_tuple_per_repeat_seq[seq_tuple.repeat_seq] = seq_tuple
    return (haps_per_repeat_seq, rep_seq_tuple_per_repeat_seq)


def read_in_asm(asm_sequences_fname):
    seq_per_hap = {}
    with open(asm_sequences_fname) as f:
        for hap_obj in SeqIO.parse(f, format="fasta"):
            if len(hap_obj.seq) < MIN_SEQ_LEN:
                continue
            hap = hap_obj.name.split('[')[0]
            assert hap not in seq_per_hap
            # trim to exact region for output
            trimmed_seq = trim_sequence(str(hap_obj.seq), config['begin_seqs'], config['end_seqs'], hap)
            if trimmed_seq is not None:
                seq_per_hap[hap] = trimmed_seq
    return seq_per_hap


def read_in_trgt(trgt_sequences_fname):
    # two sequences for each sample but unknown which hap it is
    seq_per_hap = {}
    with open(trgt_sequences_fname) as f:
        f.readline()
        for i, line in enumerate(f):
            line = line.strip()
            #print(line)
            sample, seq = line.split()
            hap = sample + ('_hap1' if (i % 2 == 0) else '_hap2')
            assert hap not in seq_per_hap
            # trim to exact region for output
            trimmed_seq = trim_sequence(seq, config['begin_seqs'], config['end_seqs'], hap)
            # **Note these haplotypes don't necessarily equate to ASM _hap1 and _hap2
            if trimmed_seq is not None:
                seq_per_hap[hap] = trimmed_seq
    return seq_per_hap


def process_sequences(
    unique_seqs, haps_per_repeat_seq, rep_seq_tuple_per_repeat_seq,
    sort_by='hap', min_haps=0, min_repeats=0, seqs_to_exclude=None, trgt_haps=None,
    noncanonical_interrs_only=False
):
    assert sort_by in ['hap', 'repeat']
    out_lines = []
    interruptions_d_overall = defaultdict(int)  # total interruptions incl multiple per hap (to print)
    interruptions_d_overall_num_haps = defaultdict(int) # num haplotypes that HAVE interruption (to print)
    repeats_hist = defaultdict(int)
    subs_used = set([config['main_repeat_unit']])
    n_have_noncanonical_interrs = 0

    repeat_seqs = set()  # does include seqs in seqs_to_exclude
    for repeat_seq in unique_seqs:
        samples_w_hap = haps_per_repeat_seq[repeat_seq]
        num_haps = len(samples_w_hap)
        if num_haps <= 0:
            #print(f'Warning: Sequence with frequency {num_haps}: {repeat_seq}') # should print only if male w two diff haplotypes
            continue
        if num_haps < min_haps:
            continue
        if trgt_haps:  # Add # TRGT haps for matching seqs output
            num_trgt_haps = 0
            if repeat_seq in trgt_haps:
                num_trgt_haps = len(trgt_haps[repeat_seq])

        (trimmed_seq, repeat_seq, seq_pre_repeat, seq_post_repeat) = rep_seq_tuple_per_repeat_seq[repeat_seq]

        interr_region = repeat_seq[repeat_seq.index(config['main_repeat_unit']):repeat_seq.rindex(config['main_repeat_unit']) + len(config['main_repeat_unit'])] # actually interrupting repeats
        repeat_seqs.add(repeat_seq)
        if (seqs_to_exclude) and (repeat_seq in seqs_to_exclude):
            continue
        out_seq_repeat_patterns = repeat_seq
        out_seq_repeat_patterns_no_tags = repeat_seq
        out_seq_repeat_patterns_maintain_num_chars = repeat_seq
        interruptions_d = defaultdict(int)
        interrs_w_ind = []
        for sub in config['substitutions']:
            expr = config['special_expr_substitutions'][sub] if (sub in config['special_expr_substitutions']) else sub 
            curr_matches = [(z.start(), sub) for z in re.finditer(expr, out_seq_repeat_patterns_maintain_num_chars) if sub != config['main_repeat_unit']]
            interrs_w_ind += curr_matches
            if curr_matches:
                subs_used.add(sub)
            out_seq_repeat_patterns = re.sub(expr,
                f'<span style="color: {config['color_per_repeat_unit'][sub]};">{REPEAT_UNIT_SHAPE}</span>',
                out_seq_repeat_patterns
            )
            out_seq_repeat_patterns_no_tags = re.sub(expr, REPEAT_UNIT_SHAPE, out_seq_repeat_patterns_no_tags)
            out_seq_repeat_patterns_maintain_num_chars = re.sub(expr, REPEAT_UNIT_SHAPE * len(sub), out_seq_repeat_patterns_maintain_num_chars)

        # interruptions processing
        for i, interr in sorted(interrs_w_ind):
            interruptions_d[interr] += 1
            interruptions_d_overall[interr] += num_haps
        interruptions_out_str = ','.join([f'{z}:{str(interruptions_d[z])}' for z in interruptions_d])

        # num haps containing interruption
        for i, interr in set(interrs_w_ind):
            interruptions_d_overall_num_haps[interr] += num_haps
        if not interrs_w_ind:
            interruptions_d_overall_num_haps[None] += num_haps

        has_noncanonical_interrs = any([z[1] not in config['interruptions_to_count'] for z in set(interrs_w_ind)])
        if has_noncanonical_interrs:
            n_have_noncanonical_interrs += 1
        if noncanonical_interrs_only and (not has_noncanonical_interrs):
            continue
        out_seq = seq_pre_repeat + out_seq_repeat_patterns + seq_post_repeat
        out_seq_no_tags = seq_pre_repeat + out_seq_repeat_patterns_no_tags + seq_post_repeat
        num_repeats = interr_region.count(config['main_repeat_unit']) + sum([interr_region.count(z) for z in config['interruptions_to_count']])
        if num_repeats < min_repeats:
            continue
        repeats_hist[num_repeats] += num_haps

        out_line = OutputLine(
            seq=f'{out_seq}',
            num_repeats=num_repeats,
            num_haps=num_haps,
            interrs_str=f'{interruptions_out_str}',
            padding_i=len(seq_pre_repeat + out_seq_repeat_patterns),
            seq_no_tags=out_seq_no_tags,
            num_trgt_haps=num_trgt_haps if (trgt_haps) else None
        )
        out_lines.append(out_line)

    unused_subs = [z for z in config['substitutions'] if z not in subs_used]
    if unused_subs:
        print('Unused encodings: ' + ', '.join(unused_subs))
    print(dict(sorted(repeats_hist.items())))

    if sort_by == 'hap':
        out_lines = sorted(out_lines, key=lambda z:(-(z[2]), z[1], z))
    else:
        out_lines = sorted(out_lines, key=lambda z:(z[1], -(z[2]), z))

    print(f'{n_have_noncanonical_interrs} alleles have noncanonical interruptions')

    return(out_lines, subs_used, repeat_seqs)


def trim_sequence(seq, begin_seqs, end_seqs, hap=None):
    repeat_begin_i = None
    for x in begin_seqs:
        matches = list(re.finditer(x, seq))
        if matches:
            m = matches[0]
            seq = seq[m.start():]  # trim beginning
            repeat_begin_i = m.end() - m.start()
            break
    if repeat_begin_i is None:
        print('Skipping- begin ', hap, seq)
        return None

    repeat_end_i = None
    for x in end_seqs:
        matches = list(re.finditer(x, seq))
        if matches:
            m = matches[0]
            seq = seq[:m.end()]  # trim end
            repeat_end_i = m.start()
            break
    if repeat_end_i is None:
        print('Skipping- end ', hap, seq)
        return None

    return TrimmedSequence(
        trimmed_seq=seq,
        repeat_seq=seq[repeat_begin_i:repeat_end_i],
        seq_pre_repeat=seq[:repeat_begin_i],
        seq_post_repeat=seq[repeat_end_i:]
    )


def write_html_divided(
    out_fname, out_lines, subs_used, divide_into=1
):
    #print(f'Total lines: {len(out_lines)}')
    out_index = list(range(0, len(out_lines), math.ceil(len(out_lines) / divide_into)))
    part_num = 1
    max_haps_overall = max([z.num_haps for z in out_lines])
    if len(out_index) > 1:
        for i in range(1, len(out_index)):
            #print(f'Writing lines {out_index[i-1]} to {out_index[i]}')
            write_html(
                out_fname.replace('.html', f'_part{part_num}.html'),
                out_lines[out_index[i-1]:out_index[i]], subs_used, write_legend=False,
                max_haps=max_haps_overall
            )
            part_num += 1
    #print(f'Writing lines {out_index[-1]} to end')
    out_fname_final = out_fname.replace('.html', f'_part{part_num}.html') if (part_num > 1) else out_fname
    write_html(out_fname_final, out_lines[out_index[-1]:], subs_used, max_haps=max_haps_overall)


def write_html(out_fname, out_lines, subs_used, write_legend=True, max_haps=0):
    if write_legend:
        out_legend = str(out_fname).replace('.html', '_legend.html')
    max_seq_len = 0
    max_trgt_haps = 0
    for l in out_lines:
        max_seq_len = max(len(re.sub(f'[^A-Z{REPEAT_UNIT_SHAPE}]', '', l.seq)), max_seq_len)
        max_haps = max(max_haps, l.num_haps)
        if l.num_trgt_haps is not None:
            max_trgt_haps = max(max_trgt_haps, l.num_trgt_haps)

    doc, tag, text = Doc().tagtext()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            with tag('title'):
                text(f'{config['gene_name']} Sequences')
            with tag('style', type='text/css'):
                text('body { font-family: "DejaVu Sans Mono", monospace; }')
                text('th, td { padding-right: 8px; text-align: left; }')
    with tag('body'):
        with tag('table'):
            with tag('thead'):
                with tag('tr'):
                    with tag('th'):
                        text('Sequence')
                    with tag('th'):
                        text('#Reps')
                    with tag('th'):
                        text('#Haps')
                    if out_lines and out_lines[0].num_trgt_haps is not None:
                        with tag('th'):
                            text('#TRGT_Haps')
                    with tag('th'):
                        text('Interruptions')
                    #with tag('th'):
                    #    text('Matches TRGT')
                    #with tag('th'):
                    #    text('TRGT Sequences (Matching if Possible)')
            with tag('tbody'):
                for out_line in out_lines:
                    with tag('tr'):
                        if out_line.num_trgt_haps is None:
                            out_line_visible = [
                                out_line.seq, out_line.num_repeats, out_line.num_haps, out_line.interrs_str
                            ]
                        else:
                            out_line_visible = [
                                out_line.seq, out_line.num_repeats, out_line.num_haps,
                                out_line.num_trgt_haps, out_line.interrs_str
                            ]
                        for i, val in enumerate(out_line_visible):
                            with tag('td'):
                                if i == 0: # Sequence
                                    val = val[:out_line.padding_i] + (
                                        '-' * (max_seq_len - len(out_line.seq_no_tags))
                                    ) + val[out_line.padding_i:]
                                    doc.asis(str(val))
                                elif i == 2 or (out_line.num_trgt_haps is not None and i == 3): # num haplotypes
                                    if max_haps >= MIN_HAPS_SHOW_BARS:
                                        with tag('div', style=f'width: {(val / max_haps) * 47}px; height: 15px; background-color: #d4b791;'):
                                            doc.asis('&nbsp;'+str(val))
                                    else:
                                        doc.asis(str(val))
                                else:
                                    doc.asis(str(val))
                        text('\n')
    with open(out_fname, 'w') as out_f:
        out_f.write(doc.getvalue())


    if write_legend:
        doc, tag, text = Doc().tagtext()
        doc.asis('<!DOCTYPE html>')
        with tag('html'):
            with tag('head'):
                with tag('title'):
                    text(f'{config['gene_name']} Legend')
                with tag('style', type='text/css'):
                    text('body { font-family: "DejaVu Sans Mono", monospace; }')
                    text('th, td { padding-right: 8px; text-align: left; }')

            with tag('body'):
                with tag('table'):
                    legend_first = [config['main_repeat_unit']] + config['interruptions_to_count']
                    for seq in legend_first + [z for z in config['color_per_repeat_unit'] if (z not in legend_first and z in subs_used)]:
                        color = config['color_per_repeat_unit'][seq]
                        with tag('tr'):
                            with tag('td'):
                                text(f'{seq}')
                            with tag('td', style=f'color: {color};'):
                                text(REPEAT_UNIT_SHAPE)
        with open(out_legend, 'w') as out_f:
            out_f.write(doc.getvalue())



if __name__ == '__main__':
    main()

