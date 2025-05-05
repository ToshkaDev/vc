#!/bin/bash

# Broad Institute's compatable files:
GENOME_LINK=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.gz
GENOME_ANNOTATION_LINK=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz

DATA_DIR=/data
GENOME_DIR=/data/genome
REFERENCE=genome.fa
ANNOTATION=genome.gtf

# [I commented out "Obtain the genome" section as the devlopment is just with the chromosome 22
# while the entire chromosome is too big. Uncomment the below if your machine
# has enough RAM (~100Gb) to work with the entire human genome] 

# Obtain the genome.
#echo "Downloading compressed genome in fasta format ..."
#wget -O .${GENOME_DIR}/${REFERENCE}.gz ${GENOME_LINK}
#echo "Unpacking genome from archive..."
#gunzip -c .${GENOME_DIR}/${REFERENCE}.gz > .${GENOME_GENOME_DIRFOLDER}/${REFERENCE}
#rm .${GENOME_DIR}/${REFERENCE}.gz

# Obtain the genome annotation:
echo "Downloading genome annotation file ..."
wget -O .${GENOME_DIR}/${ANNOTATION}.gz ${GENOME_ANNOTATION_LINK}
echo "Unpacking genome annotation from archive ..."
gunzip -c .${GENOME_DIR}/${ANNOTATION}.gz > .${GENOME_DIR}/${ANNOTATION}
rm .${GENOME_DIR}/${ANNOTATION}.gz

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
