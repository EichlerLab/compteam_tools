#!/bin/env bash

### Run this from the directory you want the Flagger jobs to run in.
# Argument is a config.sh file. Example config.sh:
#  SAMPLES_TO_RUN="samples.txt"
#  FLAGGER_PARAMS_JSON="flagger_params_w_sample_placeholder.json"
#  JSON_SAMPLE_PLACEHOLDER="insert_sample_here"  # use string not found elsewhere in json
#  FLAGGER_PATH="/net/eichler/vol26/7200/software/pipelines/flagger-0.3.3/flagger_sge.sh"
#
# Example samples.txt:
#  HG01457-hifiasm
#  HG02011-hifiasm
######

config_file="$1"
source $(pwd)/"$config_file"

######
# Launch end-to-end Flagger job for each sample
######
module load cromwell/85
mkdir -p temp
mkdir -p out_jsons
ln -sf "${FLAGGER_PATH}" $( basename "${FLAGGER_PATH}" )
while read -r sample; do
    sed "s/${JSON_SAMPLE_PLACEHOLDER}/${sample}/g" ${FLAGGER_PARAMS_JSON} > temp/flagger_params_${sample}.json
    echo "${FLAGGER_PATH} temp/flagger_params_${sample}.json out_jsons/output_${sample}.json" | qsub -V -cwd -j y -o ./log -e ./log -N run_flagger_end_to_end_${sample} -l h_rt=240:00:00 -l mfree=16G -pe serial 1 -w n -S /bin/bash
done < "${SAMPLES_TO_RUN}"