#!/bin/env bash

### For cromwell-executions directories that have final beds (not running/ failed),
# copy over final bed and remove directory.
# Run this from directory containing cromwell-executions.
# ** Note ** This will delete whole cromwell-executions directories minus final bed

mkdir -p final_beds
found_any_outfile=false
for j in $(ls cromwell-executions/*FlaggerEndToEnd*); do
    found_outfile=false
    # nohap outpath has the correct sample name in it, but prediction bed is desired output with all blocks including Hap- use this and rename
    nohap_outpath=$(find cromwell-executions/*FlaggerEndToEnd*/"$j" -path '*/output/*flagger.no_Hap.bed')
    full_outpath=$(find cromwell-executions/*FlaggerEndToEnd*/"$j" -path '*/output/*prediction.bed')
    if [ -f "$nohap_outpath" ] && [ -f "$full_outpath" ]; then
        mv "$full_outpath" final_beds/$(basename "$nohap_outpath")
        found_outfile=true
        found_any_outfile=true
    fi

    if $found_outfile; then
        echo Removing cromwell-executions/*FlaggerEndToEnd*/"$j"
        rm -rf cromwell-executions/*FlaggerEndToEnd*/"$j"
    else
        echo Skipping "$j"- Missing final bed
    fi
done

if $found_any_outfile; then
    echo "Wrote final beds to final_beds directory"
else
    echo "No output files found"
fi
