#!/bin/env bash

### For cromwell-executions directories that have final beds (not running/ failed),
# copy over final bed and remove directory.
# Run this from directory containing cromwell-executions.
# ** Note ** This will delete whole cromwell-executions directories minus final bed

mkdir -p final_beds_alt_removed

for j in $(ls cromwell-executions/*FlaggerEndToEnd*); do
    outpath=$(find cromwell-executions/*FlaggerEndToEnd*/"$j" -path '*/output/*final.bed' -o -path '*/output/*flagger.no_Hap.bed' 2>/dev/null)
    if [ -f "$outpath" ]; then
        echo Removing cromwell-executions/*FlaggerEndToEnd*/"$j"
        mv "$outpath" final_beds_alt_removed
        rm -rf cromwell-executions/*FlaggerEndToEnd*/"$j"
    else
        echo Skipping "$j"- Missing final bed
    fi
done