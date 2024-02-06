#!/bin/env bash

find -type f | sed 's/.\///' | xargs -i md5sum {} > files.md5
find -type f | sed 's/.\///' | xargs -i /bin/ls -l {} | awk '{print $5"\t"$NF}' > files.size

for file in $( cut -f 2 files.size ); do
    paste <( grep -w ${file}$ files.size ) <( grep -w ${file}$ files.md5 | awk '{print $1}' ) | awk -vOFS="\t" '{print $2,$1,$3}' 
done > manifest.txt

grep -vEw "files.size|files.md5" manifest.txt > $( basename $( pwd ) )_manifest.txt

rm files.md5
rm files.size
rm manifest.txt
