'''
Make a figure where colored rectangles for haplotypes = 100bp segments aligned more
closely to cyp 2d6/ 7, or neither/ both
'''
import pysam
import sys
from collections import defaultdict
from yattag import Doc

sam_in = sys.argv[1]  # trimmed CYP2D6-7 + flanking sample sequences divided into WINDOW_SIZE bp k-mers, aligned to fasta with CYP2D6 and CYP2D7 reference seqs
out_fname = sys.argv[2]  # out HTML

WINDOW_SIZE = 100
WINDOW_SHAPE = '▉'
WINDOW_COLORS = {
    'CYP2D6':'gold',
    'CYP2D7':'dodgerblue',
    'REP':'pink',
    'Spacer': 'tan',
    'Both':'green',
    'Neither':'lightgray',
    'Past_end':'ghostwhite'
}

max_window = 0
all_sample_haps = set()
sam_f = pysam.AlignmentFile(sam_in)
# get maximum window
for r in sam_f:
    window_start = int(r.qname.split('_')[-3])
    max_window = max(max_window, window_start)
    all_sample_haps.add('_'.join(r.query_name.split('_')[:2]))
all_windows = range(0, max_window, WINDOW_SIZE)

mapping_d = {}
for sample_hap in all_sample_haps:
    if sample_hap not in mapping_d:
        mapping_d[sample_hap] = {}

sam_f = pysam.AlignmentFile(sam_in)
for r in sam_f:
    sample_hap = '_'.join(r.query_name.split('_')[:2])
    mapq = int(r.mapq)
    flag = r.flag
    window_start = int(r.qname.split('_')[-3])
    if window_start not in mapping_d[sample_hap]:
        mapping_d[sample_hap][window_start] = set()
    if r.flag & 0x4 == 0: # is mapped
        gene = r.reference_name.split('_')[0]
        mapping_d[sample_hap][window_start].add((gene, mapq))

renamed_samples = {z:str(i + 1) for i, z in enumerate(set([z.split('_')[0] for z in all_sample_haps]))}
renamed_haps_ordered = {}
for x in renamed_samples:
    if f'{x}_hap1' in all_sample_haps:
        renamed_haps_ordered[f'{x}_hap1'] = f'{renamed_samples[x]}_hap1' # f'{x}_hap1'
    if f'{x}_hap2' in all_sample_haps:
        renamed_haps_ordered[f'{x}_hap2'] = f'{renamed_samples[x]}_hap2' # f'{x}_hap2' 

out_lines = []
for sample_hap in renamed_haps_ordered:
    hap_d = mapping_d[sample_hap]
    out_str = ''
    for window in all_windows:
        best_mapping_genes = []
        if window in hap_d:
            if hap_d[window]:
                max_gene_score = max([z[1] for z in hap_d[window]])
                best_mapping_genes = [z[0] for z in hap_d[window] if z[1] == max_gene_score]
                color_key = None
                if 'CYP2D6' in best_mapping_genes and 'CYP2D7' in best_mapping_genes:
                    color_key = 'Both'
                elif 'CYP2D6' in best_mapping_genes:
                    color_key = 'CYP2D6'
                elif 'CYP2D7' in best_mapping_genes:
                    color_key = 'CYP2D7'
                elif 'REP6' in best_mapping_genes:
                    color_key = 'REP'
                elif 'Spacer' in best_mapping_genes:
                    color_key = 'Spacer'
                else:
                    color_key = 'Neither'
            else:
                color_key = 'Neither'
        else:
            color_key = 'Past_end'
        out_str += f'<span style="color: {WINDOW_COLORS[color_key]};">{WINDOW_SHAPE}</span>'
    out_lines.append((renamed_haps_ordered[sample_hap], out_str))
# sort by most mapped blocks
out_lines = sorted(out_lines, key=lambda x:x[1].count(WINDOW_COLORS['CYP2D6']) + x[1].count(WINDOW_COLORS['CYP2D7']) + x[1].count(WINDOW_COLORS['Both']), reverse=True)

doc, tag, text = Doc().tagtext()
doc.asis('<!DOCTYPE html>')
with tag('html'):
    with tag('head'):
        with tag('title'):
            text('CYP2D6/7 Mappings')
        with tag('style', type='text/css'):
            text('body { font-family: "DejaVu Sans Mono", monospace; }')
            text('th, td { padding-right: 8px; text-align: left; }')
    with tag('body'):
        with tag('table'):
            with tag('thead'):
                with tag('tr'):
                    with tag('th'):
                        text('Renamed_Sample')
                    with tag('th'):
                        text('Better-Mapping_Gene')
            with tag('tbody'):
                for out_line in out_lines:
                    with tag('tr'):
                        with tag('td', style='white-space: nowrap;'):
                            doc.asis(f'{out_line[0]}')
                        with tag('td'):
                            doc.asis(str(out_line[1]))
                    text('\n')
    doc.stag('hr')
    doc.asis('<b>Assembly Position Maps Better To:</b>')
    with tag('table'):
        for gene in WINDOW_COLORS:
            color = WINDOW_COLORS[gene]
            with tag('tr'):
                with tag('td'):
                    text(f'{gene}')
                with tag('td', style=f'color: {color};'):
                    text(WINDOW_SHAPE)

    with open(out_fname, 'w') as out_f:
        out_f.write(doc.getvalue())

