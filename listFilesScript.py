import sys
import os
import pandas as pd
import hashlib

#root = "/net/eichler/vol28/home/nidhi12k/nobackups/melanie_test/"
root=sys.argv[1]

def get_file_size_gb(file_path):
    file_size_bytes = os.path.getsize(file_path)
    file_size_gb = file_size_bytes / (1024 ** 3)
    return file_size_bytes

def calculate_md5(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
        return hash_md5.hexdigest()

def get_md5(file_path):
    with open(file_path, "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)
        return file_hash.hexdigest()

df = pd.DataFrame()
list_filenames = []
list_topdir = []
list_size = []
list_md5sum = []
list_file_path = []

for path, subdirs, files in os.walk(root):
    for name in files:
        if "hifi_reads" not in name:
            continue
        elif name.endswith('.bam'):
            file_path = os.path.join(path, name)
            rel_path = os.path.relpath(file_path, root)
            top_dir = rel_path.split(os.sep)[0]
            file_size = get_file_size_gb(file_path) 
            list_filenames.append(name)
            list_topdir.append(top_dir)
            list_size.append(file_size)
            list_file_path.append(file_path)

df = pd.DataFrame({
	'filename': list_filenames,
	'top_level_dir': list_topdir,
	'size': list_size,
	'file_path': list_file_path
})

rows_to_remove = df[df['filename'].str.contains("unassigned")].index
df = df.drop(rows_to_remove)
rows_to_remove = df[df['filename'].str.contains("fail_reads")].index
df = df.drop(rows_to_remove)

#df["md5sum"] = df["file_path"].apply(get_md5)
#df = df[['filename', 'top_level_dir', 'size', 'md5sum', 'file_path']]
df = df[['filename', 'top_level_dir', 'size', 'file_path']]

df.to_csv("file_list.tsv", sep="\t", index=None)
