#!/bin/bash

set -x
set -e

## Configurable file paths
GUIDEDIR="../../library_v1/"
BARCODINGDIR="../../barcode_v1_181005/"

WORKDIR="./work"
FIGUREDIR="./figures"

FASTQDIR="/mnt/ingolialab/FastQ/190425_50SR_HS4K2A/Meacham"

BARCODE_ASSIGN_DIR="../../barcode_v1_181005/assignment/barcode-assign-0.1"
BARCODE_ASSIGN_PATH="${BARCODE_ASSIGN_DIR}/target/debug/"

BC_COUNT="${BARCODE_ASSIGN_PATH}/bc-count"
BC_TABULATE="${BARCODE_ASSIGN_PATH}/bc-tabulate"

CONSTANT="GCATGCGTGAAGTGGCGCGCCTGATA"

if [[ ! -x "${BC_COUNT}" || ! -x "${BC_TABULATE}" ]];
then
    echo "bc-count and bc-tabulate needed"
    echo "build the barcode-assign tools in \"${BARCODE_ASSIGN_DIR}\""
    exit 1
fi

mkdir -p "${WORKDIR}"
mkdir -p "${FIGUREDIR}"

COUNT_LIST=""

for FASTQFILE in `cut -f1 samples.txt`
do
    SAMPLE=`grep -F ${FASTQFILE} samples.txt | cut -f2`

    BASE="${WORKDIR}/${SAMPLE}"
    COUNT="${BASE}-nbhd-count.txt"
    
    if [[ ! -e "${COUNT}" ]];
    then
        cutadapt -a "${CONSTANT}" \
	       --minimum-length 10 --discard-untrimmed \
	       "${FASTQDIR}/${FASTQFILE}" 2>"${WORKDIR}/${SAMPLE}-cutadapt.txt" \
	  | "${BC_COUNT}"  -f - -o "${BASE}-nbhd-count.txt" -n "${BASE}" &
        sleep 1
    else
        echo "${COUNT} exists, skipping ${SAMPLE}..."
    fi

    COUNT_LIST="${COUNT_LIST} ${COUNT}"
done

wait

if [[ ! -e "${WORKDIR}/nizm005.txt" ]];
then
    "${BC_TABULATE}" --minsamples 2 -o "${WORKDIR}/nizm005.txt" ${COUNT_LIST}
fi

export GUIDEDIR
export BARCODINGDIR
export WORKDIR
export FIGUREDIR

if [[ ! -e "${WORKDIR}/nizm005-assigned.txt" ]];
then  
    R --no-save < analyze_barcodes.R
else
    echo "\"${WORKDIR}/nizm005-assigned.txt\" exists, skipping analyze_barcodes.R"
fi

if [[ ! -e "${FIGUREDIR}/nizm005-guide-deseq.csv" ]];
then
    R --no-save < deseq_barcodes.R
else
    echo "\"${FIGUREDIR}/nizm005-guide-deseq.csv\" exists, skipping deseq_barcodes.R"
fi

if [[ ! -e "${WORKDIR}/nizm005-yorf-deseq.csv" ]];
then
    R --no-save < analyze_fitness.R
else
    echo "\"${WORKDIR}/nizm005-yorf-deseq.csv\" exists, skipping analyze_fitness.R"
fi

if [[ ! -e "${FIGUREDIR}/guide-table-nizm005-deseq.csv" ]];
then
    R --no-save < tabulate_guides.R
else
    echo "\"${FIGUREDIR}/guide-table-nizm005-deseq.csv\" exists, skipping tabulate_guides.R"
fi

R --no-save < analyze_guides.R
