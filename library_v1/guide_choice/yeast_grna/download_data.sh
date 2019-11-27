#!/bin/bash

set -e
set -x

if [[ -z "$1" ]];
then
    echo "Usage: $0 <DATADIR>"
    exit 1
fi

export DATADIR="$1"

if [[ ! -d "${DATADIR}" ]];
then
    echo "Data directory \"${DATADIR}\" does not exist"
    exit 1
fi

GENOMEFA="${DATADIR}/sacCer3.fa"
GENOMEBT2="${DATADIR}/sacCer3"

## Yeast genome Fasta files
## One file per chromosome in gzipped tar
## chrom.sizes used for a table of chromosome names
CHROMGZ="${DATADIR}/sacCer3.tar.gz"
GENOMEDIR="${DATADIR}/sacCer3"
CHROMSIZES="${GENOMEDIR}/chrom.sizes.txt"

## Build a single genomic Fasta file, one entry per chromosome
if [[ ! -e "${GENOMEFA}" ]];
then
    rm -f "${GENOMEFA}"

    if [[ ! -e "${CHROMGZ}" ]];
    then
        curl 'http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz' --output "${CHROMGZ}"
    else
        echo "${CHROMGZ} exists, not downloading..."
    fi
    
    mkdir -p "${GENOMEDIR}"
    tar -xv -C "${GENOMEDIR}" -f "${CHROMGZ}"
    
    if [[ ! -e "${CHROMSIZES}" ]];
    then
        curl 'http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes' --output "${CHROMSIZES}"
    else
        echo "${CHROMSIZES} exists, not downloading..."
    fi
    
    CHROMS=`cut -f1 "${CHROMSIZES}"`

    for CHROM in ${CHROMS}
    do
        CHROMFA="${GENOMEDIR}/${CHROM}.fa"
        if [[ ! -e "${CHROMFA}" ]];
        then
	  echo "chr.fa file ${CHROMFA} for ${CHROM} does not exist"
	  exit 1
        fi

        cat "${CHROMFA}" >> "${GENOMEFA}"

    done
else
    echo "${GENOMEFA} exists, not compiling"
fi

## Bowtie2 index of genome Fasta file
if [[ "${GENOMEFA}" -nt "${GENOMEBT2}.1.bt2" ]];
then
    bowtie2-build "${GENOMEFA}" "${GENOMEBT2}"
else
    echo "Skipping bowtie2 index, ${GENOMEBT2}.1.bt2 exists newer than ${GENOMEFA}"
fi

## BED-format annotations
## UCSC sgdGene table is close to BED
GENOMEBED="${DATADIR}/sacCer3.bed"
SGDGENE="${DATADIR}/sgdGene.txt"
SGDGENEGZ="${DATADIR}/sgdGene.txt.gz"

if [[ ! -e "${GENOMEBED}" ]];
then
    if [[ ! -e "${SGDGENE}" ]];
    then
        curl "http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/database/sgdGene.txt.gz" -o "${SGDGENEGZ}"
        gunzip "${SGDGENEGZ}"
    fi

    awk -F$'\t' '{printf("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0\t%d\t%s\t%s\n", $3, $5, $6, $2, $4, $7, $8, $9, $10, $11)}' < "${SGDGENE}" > "${GENOMEBED}"
fi

## SGD features indicate which ORFs are "dubious"
DUBIOUS="${DATADIR}/dubious.txt"
SGDFEATURES="${DATADIR}/SGD_features.tab"
if [[ ! -e "${DUBIOUS}" ]];
then
    if [[ ! -e "${SGDFEATURES}" ]];
    then
        curl 'https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab' -o "${SGDFEATURES}"
    fi

    awk -F$'\t' '{ if ($3 == "Dubious") { print $4 } }' < "${SGDFEATURES}" | sort > "${DUBIOUS}"
else
    echo "Dubious gene list ${DUBIOUS} exists..."
fi

## Extract 5' end from TIF

## Pelechano, V., Wei, W., Steinmetz, L.M., 2013. Extensive
## transcriptional heterogeneity revealed by isoform profiling. Nature
## 497, 127–131. doi:10.1038/nature12121

ZIP2="${DATADIR}/Supplementary-Data-2.zip"
ANNO_MTIF="${DATADIR}/S2_tcd_mTIFAnno2.txt"

# We clustered the transcripts with 5′ and 3′ end sites co-occurring
# within 5 bp (Supplementary Fig. 2a). Specifically, we defined TSS and
# TTS clusters separately. Each cluster was defined by both a window and
# a mTSS/TTS (the most abundant within that window). Clusters of
# TSS/TTSs were assigned iteratively in decreasing order of
# expression...After this assignment process, only clusters defined by
# mTSS/TTSs with at least 3 supporting reads were considered. mTIFs were
# defined as connections between mTSS and mTTS supported by at least 2
# reads connecting the associated clusters. All TIFs that shared a given
# TSS/TTS cluster were assigned to the corresponding mTIF cluster.

if [[ ! -e "${ANNO_MTIF}" ]];
then    
    if [[ ! -e "${ZIP2}" ]];
    then
        curl 'http://steinmetzlab.embl.de/TIFSeq/data/Supplemental/Supplementary%20Data%202.zip' --output "${ZIP2}"
    fi

    unzip "${ZIP2}" -d `dirname ${ANNO_MTIF}` `basename ${ANNO_MTIF}`
fi

## ATAC-seq data from GEO
## Schep AN, Buenrostro JD, Denny SK, Schwartz K, Sherlock G, Greenleaf WJ.
## Structured nucleosome fingerprints enable high-resolution mapping of chromatin
## architecture within regulatory regions. Genome Res. 2015 Nov;25(11):1757-70. doi:
## 10.1101/gr.192294.115.
GZFILE="${DATADIR}/GSE66386_all_occ.bw.tar.gz"
GSEDIR="${DATADIR}/GSE66386"
OCCFILES="A.occ.wig B.occ.wig"

if [[ ! -e "${GZFILE}" ]];
then
    set +e
    curl 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE66386&format=file&file=GSE66386%5Fall%5Focc%2Ebw%2Etar%2Egz' --output "${GZFILE}"
    set -e
    if [[ ! -e "${GZFILE}" ]];
    then
	echo "Failed to download ${GZFILE}"
	exit 1
    fi
else
    echo "${GZFILE} exists, not downloading..."
fi

mkdir -p "${GSEDIR}"

for WIG in ${OCCFILES}
do
    WIGPATH="${GSEDIR}/${WIG}"

    if [[ ! -e "${WIGPATH}" ]];
    then
        BASE=`basename ${WIG} .wig`
        BWPATH="${GSEDIR}/${BASE}.bw"
        if [[ ! -e "${BWPATH}" ]];
        then
	  tar -xv -C "${GSEDIR}" -f "${GZFILE}" "${BASE}.bw"
        fi
        bigWigToWig "${GSEDIR}/${BASE}.bw" "${GSEDIR}/${WIG}"
    else
        echo "${WIGPATH} exists, not converting..."
    fi
done

