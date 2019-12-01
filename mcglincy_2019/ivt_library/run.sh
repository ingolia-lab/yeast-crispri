#!/bin/bash

set -x
set -e

## Configurable file paths
DATADIR="./work"
FASTQDIR="/mnt/ingolialab/FastQ/181129_50SR_HS4KB/Meacham/"

BARCODE_ASSIGN_DIR="../../barcode_v1_181005/assignment/barcode-assign-0.1"
BARCODE_ASSIGN_PATH="${BARCODE_ASSIGN_DIR}/target/debug/"

BC_COUNT="${BARCODE_ASSIGN_PATH}/bc-count"
BC_TABULATE="${BARCODE_ASSIGN_PATH}/bc-tabulate"

# 3' adapter sequence for barcode
CONSTANT="GCATGCGTGAAGTGGCGCGCCTGATA"

if [[ ! -x "${BC_COUNT}" || ! -x "${BC_TABULATE}" ]];
then
    echo "bc-count and bc-tabulate needed"
    echo "build the barcode-assign tools in \"${BARCODE_ASSIGN_DIR}\""
    exit 1
fi

mkdir -p "${DATADIR}"

COUNT_LIST=""

for FASTQFILE in `cut -f1 samples.txt`
do
    SAMPLE=`grep -F ${FASTQFILE} samples.txt | cut -f2`

    BASE="${DATADIR}/${SAMPLE}"
    COUNT="${BASE}-nbhd-count.txt"
    
    if [[ ! -e "${COUNT}" ]];
    then
        cutadapt -a "${CONSTANT}" \
         	       --minimum-length 10 --discard-untrimmed \
        	       "${FASTQDIR}/${FASTQFILE}" 2>"${DATADIR}/${SAMPLE}-cutadapt.txt" \
         	  | "${BC_COUNT}"  -f - -o "${BASE}-nbhd-count.txt" -n "${BASE}" &
        sleep 1
    else
        echo "${COUNT} exists, skipping ${SAMPLE}..."
    fi

    COUNT_LIST="${COUNT_LIST} ${COUNT}"
done

wait

"${BC_TABULATE}" --minsamples 2 -o "${DATADIR}/ivt_all.txt" ${COUNT_LIST}

export DATADIR
R --no-save < filter_ivt.R
R --no-save < analyze_ivt.R

