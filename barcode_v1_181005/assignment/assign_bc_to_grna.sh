#!/bin/bash

set -x
set -e

## Configurable file paths

## Reference file for library oligos
TARGET_OLIGOS=../../library_v1/sequence-oligos.fa

## Paired-end sequencing files
ASSIGN_R1="/mnt/ingolialab/FastQ/181219_75PE_MS2/Meacham/NIZM002_S1_L001_R1_001.fastq.gz"
ASSIGN_R2="/mnt/ingolialab/FastQ/181219_75PE_MS2/Meacham/NIZM002_S1_L001_R2_001.fastq.gz"

## Working directory
DATADIR="./work"
OUTDIR="./"

## barcode-assign tools
## source directory needed only if binaries not present in path
BARCODE_ASSIGN_DIR="./barcode-assign-0.1"
BARCODE_ASSIGN_PATH="${BARCODE_ASSIGN_DIR}/target/debug/"

if [[ ! -x "${BARCODE_ASSIGN_PATH}/bc-seqs" ]];
then
    echo "Cannot find barcode-assign tools in \"${BARCODE_ASSIGN_PATH}\""
    if [[ ! -e "${BARCODE_ASSIGN_DIR}/Cargo.toml" ]];
    then
        echo "Cannot find barcode-assign source in \"${BARCODE_ASSIGN_DIR}\""
        exit 1
    fi
    ## N.B. 1.39 breaks rust-htslib -- waiting for fix
    (cd "${BARCODE_ASSIGN_DIR}"; rustup run 1.38.0 cargo build)
    if [[ ! -x "${BARCODE_ASSIGN_PATH}/bc-seqs" ]];
    then
        echo "Failed to build barcode-assign tools from source"
        exit 1
    fi
fi

## Adapter sequences for barcode (3' adapter) and guide (both sides)
BCTRIM=GCATGCGTGAAGTGGCGCGCCTGATAGCTCGTTTAAACTG
GTRIM5=CGAAAC
GTRIM3=AAGTTAAAAT

mkdir -p "${DATADIR}"

# Trim barcode adapter
if [[ ! -e "${DATADIR}/bctrim_R1.fq" ]]
then
    cutadapt -a "${BCTRIM}" --pair-filter=both \
	   --untrimmed-output="${DATADIR}/no_bctrim_R1.fq" \
	   --untrimmed-paired-output="${DATADIR}/no_bctrim_R2.fq" \
	   --interleaved \
	   "${ASSIGN_R1}" "${ASSIGN_R2}" \
	   2> "${DATADIR}/bctrim-adapter.txt" \
        | cutadapt --interleaved \
	         --pair-filter=any --minimum-length=12 \
	         --too-short-output="${DATADIR}/short_bctrim_R1.fq" \
	         --too-short-paired-output="${DATADIR}/short_bctrim_R2.fq" \
	         -o "${DATADIR}/bctrim_R1.fq" -p "${DATADIR}/bctrim_R2.fq" \
	         - \
	         > "${DATADIR}/bctrim-length.txt"  2>&1
fi

# Use barcode (R1) sequence as name for guide (R2) sequence
if [[ ! -e "${DATADIR}/grna_barcode_barcoded.fq" ]]
then
    "${BARCODE_ASSIGN_PATH}/bc-seqs" \
        --barcodes "${DATADIR}/bctrim_R1.fq" \
        --sequences "${DATADIR}/bctrim_R2.fq" \
        --outbase "${DATADIR}/grna_barcode" \
        --neighborhood "${DATADIR}/grna_barcode"
fi

# Trim barcoded guide sequences
if [[ ! -e "${DATADIR}/grna_barcode_trimmed.fq" ]]
then
    cutadapt -a "${GTRIM5}...${GTRIM3}" \
	   --untrimmed-output="${DATADIR}/no_grnatrim.fq" \
	   --minimum-length=20 \
	   --too-short-output="${DATADIR}/short_grnatrim.fq" \
	   "${DATADIR}/grna_barcode_barcoded.fq" \
	   -o "${DATADIR}/grna_barcode_trimmed.fq" \
	   > "${DATADIR}/grnatrim-adapter.txt" 2>&1
fi

# Generate bowtie index of guide sequences
if [[ ! -e "${DATADIR}/target-oligos.1.bt2" ]]
then
    cp "${TARGET_OLIGOS}" "${DATADIR}/target-oligos.fa"
    echo ">Empty" >> "${DATADIR}/target-oligos.fa"
    echo "tctgggagctgcgattggGATCCTCGAAGCTttttagagctagaaatagc" >> "${DATADIR}/target-oligos.fa"
    bowtie2-build "${DATADIR}/target-oligos.fa" "${DATADIR}/target-oligos"
fi

# Align guides against bowtie index and sort output by name
if [[ ! -e "${DATADIR}/grna_barcode_sorted.bam" ]]
then
    bowtie2 -p38 \
	  -x "${DATADIR}/target-oligos" \
	  -U "${DATADIR}/grna_barcode_trimmed.fq" \
	  2> "${DATADIR}/grna-bowtie.txt" \
        | samtools view -b - 2> "${DATADIR}/grna-bowtie-samtools.txt" \
        | samtools sort -m 4G -n -o "${DATADIR}/grna_barcode_sorted.bam" - \
	         2> "${DATADIR}/grna-samtools-sort.txt"
fi

mkdir -p "${OUTDIR}"

# Determine barcode-to-guide assignments
if [[ ! -e "${OUTDIR}/grna-assign-barcode-fates.txt" ]]
then
    "${BARCODE_ASSIGN_PATH}/bc-grna" \
        --align-start 0 \
        --bam-by-name "${DATADIR}/grna_barcode_sorted.bam" \
        --outbase "${OUTDIR}/grna-assign-"
fi
