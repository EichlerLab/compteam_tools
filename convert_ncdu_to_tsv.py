import json
import sys
import os
import bz2

ncdu_fname = sys.argv[1]
out_tsv = sys.argv[2]

def main():
    with bz2.open(ncdu_fname, 'rt', errors='ignore') as f:
        j = json.load(f)
    with open(out_tsv, 'w') as out_f:
        out_f.write('Filepath\tasize\tdsize\tinode\n')
        process_entries(out_f, j[3], '')

def process_entries(out_f, entries, base_path=''):
    curr_path = os.path.join(base_path, entries[0]['name'])
    if len(entries) == 1:
        x = entries[0]
        out_f.write('\t'.join((
            curr_path, str(x.get('asize', '')), str(x.get('dsize', '')), str(x.get('ino', ''))
        )) + '\n')
        return
    for x in entries[1:]:
        if isinstance(x, list):
            process_entries(out_f, x, curr_path)
        else:
            out_f.write('\t'.join((
                os.path.join(curr_path, x['name']), str(x.get('asize', '')), str(x.get('dsize', '')), str(x.get('ino', ''))
            )) + '\n')


if __name__ == '__main__':
    main()
