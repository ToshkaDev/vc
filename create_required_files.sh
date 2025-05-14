#!/bin/bash

# Broad Institute's compatable files:
GENOME_LINK=https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
GENOME_ANNOTATION_LINK=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz

DATA_DIR=/data
GENOME_DIR=/data/genome
REFERENCE=genome.fa
ANNOTATION=genome.gtf

# [I commented out "Obtain the genome" section as the devlopment is just with the chromosome 22
# while the entire chromosome is too big. Uncomment the below if your machine
# has enough RAM (~100Gb) to work with the entire human genome] 

# Obtain the genome.
if [ ! -f .${GENOME_DIR}/${REFERENCE} ]; then 
    echo "Downloading genome in fasta format ..."
    wget -O .${GENOME_DIR}/${REFERENCE} ${GENOME_LINK}
fi

# Obtain the genome annotation:
if [ ! -f .${GENOME_DIR}/${ANNOTATION} ]; then 
    echo "Downloading genome annotation file ..."
    wget -O .${GENOME_DIR}/${ANNOTATION}.gz ${GENOME_ANNOTATION_LINK}
    echo "Unpacking genome annotation from archive ..."
    gunzip -c .${GENOME_DIR}/${ANNOTATION}.gz > .${GENOME_DIR}/${ANNOTATION}
    rm .${GENOME_DIR}/${ANNOTATION}.gz
fi

#Extracting a particular chromosome:
# docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464 bashs -c "samtools \
#     faidx Homo_sapiens_assembly38.fasta chr22 > chr22.fa"

################

# 1) Create a simple fasta index file (.fai):
echo "Creating simple fasta index file ..."
docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464 samtools \
    faidx ${REFERENCE}

# 2) Creaet a sequence dictionary used by GATK and Picard:
echo "Creating sequence dictionary ..."
docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867 gatk CreateSequenceDictionary \
    -R ${REFERENCE} \
    -O ${REFERENCE%.*}.dict

# 3) Obtaining BED file from GTF/GFF (e.g., from Ensembl/NCBI).
# BED file describes genomic intervals (e.g., genes, exons, regions of interest)
echo "Creating BED file (describes genomic intervals) ..."
docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/bedops:2.4.41--0451d22c61ea1547 bash -c "gff2bed \
< ${ANNOTATION} > ${ANNOTATION%.*}.bed"

################
# Download the RNA edit sites file for filtering called variants

RNA_EDIT_SITES_DIR=.${DATA_DIR}/rna_edit_sites
RNA_EDIT_SITES_LINK=http://srv00.recas.ba.infn.it/webshare/ATLAS/download/TABLE1_hg38_v2.txt.gz
RNA_EDIT_SITES_FILE=rna_edit_sites.gz

mkdir -p ${RNA_EDIT_SITES_DIR}
if [ ! -f ${RNA_EDIT_SITES_DIR}/${RNA_EDIT_SITES_FILE%.*}.txt ]; then
    echo "Downloading RNA edit sites file ..."
    wget -O ${RNA_EDIT_SITES_DIR}/${RNA_EDIT_SITES_FILE} ${RNA_EDIT_SITES_LINK}

    echo "Unpacking gz compressed file ... "
    gunzip  -c ${RNA_EDIT_SITES_DIR}/${RNA_EDIT_SITES_FILE} > ${RNA_EDIT_SITES_DIR}/${RNA_EDIT_SITES_FILE%.*}.txt
    rm ${RNA_EDIT_SITES_DIR}/${RNA_EDIT_SITES_FILE}
fi

################
# Download the LCR bed file for human genome for low complexity regions (LCRs) filtering

LCR_DIR=.${DATA_DIR}/lcr
LCR_LINK=https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true
LCR_FILE=lcr_with_chr.bed.gz

mkdir -p ${LCR_DIR}

if [ ! -f ${LCR_DIR}/${LCR_FILE%.*} ]; then
    echo "Downloading low complexity regions (LCR) file for a genome ..."
    curl -o ${LCR_DIR}/${LCR_FILE} -L https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true

    echo "Unpacking LCR gz compressed file ... "
    gunzip -c ${LCR_DIR}/${LCR_FILE} > ${LCR_DIR}/${LCR_FILE%.*}
    rm ${LCR_DIR}/${LCR_FILE}
fi