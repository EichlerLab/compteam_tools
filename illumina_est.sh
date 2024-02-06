#!/bin/env bash


OUT_README=$( readlink -f ${1} | xargs -i dirname {} )/COV_README

if [[ -f ${OUT_README} ]]; then
    echo -e "FOUND ${COV_README}. Not recomputing"
else
    echo -e "Estimating coverage for ${@} and leaving README in ${OUT_README}" 
    cov=$( zcat $@ | wc -l | awk '{print $1/4*151/3100000000}' )
    echo -e "Estimated coverage for ${@}: ${cov}" > ${OUT_README} 
    echo -e "Estimated coverage for ${@}: ${cov}" 
fi

