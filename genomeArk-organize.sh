#!/bin/env bash


raw_data_format="pod5"
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


####################
# PacBio HiFi data #
####################

cd ${LRA}/sharing

mkdir -p genomeark-upload

cd genomeark-upload/
mkdir -p ${GARK_NAME}/genomic_data/{pacbio_hifi,ont}

pushd ${GARK_NAME}/genomic_data/pacbio_hifi

/bin/ls /net/eichler/vol28/projects/long_read_archive/nobackups/nhp/${LRA_NAME}/raw_data/PacBio_HiFi/*bam | grep -v "\.reads\." | xargs -i cp -l {} ./
/bin/ls /net/eichler/vol28/projects/long_read_archive/nobackups/nhp/${LRA_NAME}/raw_data/PacBio_HiFi/*bam.pbi | grep -v "\.reads\." | xargs -i cp -l {} ./

/bin/ls /net/eichler/vol28/projects/long_read_archive/nobackups/nhp/${LRA_NAME}/raw_data/PacBio_HiFi/*fastq.gz | grep -v "\.reads\." | xargs -i cp -l {} ./

/bin/ls /net/eichler/vol28/projects/long_read_archive/nobackups/nhp/${LRA_NAME}/raw_data/PacBio_HiFi/subread/*subreads*bam | xargs -i cp -l {} ./
/bin/ls /net/eichler/vol28/projects/long_read_archive/nobackups/nhp/${LRA_NAME}/raw_data/PacBio_HiFi/subread/*subreads*bam.pbi | xargs -i cp -l {} ./

for file in $( /bin/ls *bam* ); do 
    paste <( grep ${file}$ <( find /net/eichler/vol28/projects/long_read_archive/nobackups/nhp/${LRA_NAME}/raw_data/PacBio_HiFi/copy_record/ -type f | xargs -i zcat {} | awk -vOFS="\t"  '{ print $12,$7 }' ) | cut -f 1 ) <( echo ${file} ) | awk '{if (NF == 2) print }'; 
done > md5.pb

for file in $( /bin/ls | grep -v ^md5 ); do
    grep -q ${file}$ md5.pb
    if [[ $? != 0 ]]; then
        md5sum ${file}
    fi
done > md5.new

cat md5.pb md5.new | awk -vOFS="  " '{print $1,$2}' > files.md5; rm md5.new md5.pb

pushd

############
# ONT data #
############

pushd ${GARK_NAME}/genomic_data/ont

find /net/eichler/vol28/projects/long_read_archive/nobackups/nhp/${LRA_NAME}/raw_data/nanopore/ -type f | grep -E "fastq.gz$|bam$" | grep pass | grep dorado | xargs -i cp -l {} ./


for dir in $( /bin/ls -d /net/eichler/vol28/projects/long_read_archive/nobackups/nhp/${LRA_NAME}/raw_data/nanopore/*/fast5/*/${raw_data_format} ); do
    runname=$( echo ${dir} | awk -F "/" '{print $(NF-1)}' ); 
    mkdir -p ${raw_data_format}/${runname};
    find ${dir} -type f | grep ${raw_data_format}$ | xargs -i cp -l {} ${raw_data_format}/${runname}
done

md5sum *.bam *.fastq.gz > files.md5

cd ${raw_data_format}
for dir in $( /bin/ls ); do
    cd ${dir}
    echo ${dir}
    md5sum *.${raw_data_format} > files.md5
    cd ../
done
cd ../

pushd
