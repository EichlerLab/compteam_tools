#!/bin/env bash

LRA="/net/eichler/vol28/projects/long_read_archive/nobackups"


while getopts g:l:f: flag
do
    case "${flag}" in
        g) GARK_NAME=${OPTARG};;
        l) LRA_NAME=${OPTARG};;
        f) raw_data_format=${OPTARG};;
    esac
done

if [[ ( $( echo ${GARK_NAME} ) == "" ) || ( $( echo ${LRA_NAME} ) == "" )  ]]; then
    echo "-g and -l are required"
    exit 1
fi

cd /net/eichler/vol28/projects/long_read_archive/nobackups/sharing/genomeark-upload/${GARK_NAME}

for DATATYPE_DIR in $( /bin/ls -d genomic_data/{pacbio_hifi,ont,arima} transcriptomic_data/brain/pacbio )
do
    nano ${DATATYPE_DIR}/README
done

UUID=$( uuidgen )

NAME=$(date +%Y%m%d)_${GARK_NAME}_ONT_HIFI_${LRA_NAME}

alias aws_gark='aws --profile genomeark'
S3GARK=s3://genomeark

aws_gark s3 sync ./ $S3GARK/incoming/${UUID}--${NAME}

