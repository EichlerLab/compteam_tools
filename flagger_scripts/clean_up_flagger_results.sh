#!/bin/env bash

### For cromwell-executions directories that have final beds (not running/ failed),
# copy over final bed and remove directory.
# Run this from directory containing cromwell-executions.
# ** Note ** This will delete whole cromwell-executions directories minus final bed

mkdir -p final_outputs
found_any_outfile=false
for j in $(ls cromwell-executions/*FlaggerEndToEnd*); do
    found_outfile=false
    # nohap outpath has the correct sample name in it, but prediction bed is desired output with all blocks including Hap- use this and rename
    nohap_outpath=$(find cromwell-executions/*FlaggerEndToEnd*/"$j" -path '*/output/*flagger.no_Hap.bed' 2>/dev/null)
    sample=$(basename "$nohap_outpath" | sed 's/\.hmm_flagger\.no_Hap\.bed//')
    full_outpath=$(find cromwell-executions/*FlaggerEndToEnd*/"$j" -path '*/output/*prediction.bed' 2>/dev/null)
    if [ -f "$nohap_outpath" ] && [ -f "$full_outpath" ]; then
        echo Sample "$sample":
        echo Copying final bed to final_outputs/"$sample"_$(basename "$full_outpath")
        cp "$full_outpath" final_outputs/"$sample"_$(basename "$full_outpath")
        found_outfile=true
        found_any_outfile=true
    else
        v033_outpath=$(find cromwell-executions/*FlaggerEndToEnd*/"$j" -path '*/output/*final.bed')
        if [ -f "$v033_outpath" ]; then
            # already has sample name
            echo Copying final bed to final_outputs/$(basename "$v033_outpath")
            cp "$v033_outpath" final_outputs
            found_outfile=true
            found_any_outfile=true
        fi
    fi

    if $found_outfile; then
        summary_outpath=$(find cromwell-executions/*FlaggerEndToEnd*/"$j" -name 'prediction_summary_final.tsv' 2>/dev/null)
        if [ -f "$summary_outpath" ]; then
            echo Copying summary TSV to final_outputs/"$sample"_$(basename "$summary_outpath")
            cp "$summary_outpath" final_outputs/"$sample"_$(basename "$summary_outpath")
        fi 

        echo Removing task directory cromwell-executions/*FlaggerEndToEnd*/"$j"
        rm -rf cromwell-executions/*FlaggerEndToEnd*/"$j"
        echo
    else
        echo Skipping "$j"- Missing final bed
        echo
    fi
done

if $found_any_outfile; then
    echo "Wrote final beds/ summary files to final_outputs directory"
else
    echo "No output files found"
fi
