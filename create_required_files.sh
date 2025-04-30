#!/bin/bash

GENOME_LINK=https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
GENOME_ANNOTATION_LINK=https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

DATA_DIR=/data
GENOME_DIR=/data/genome
REFERENCE=genome.fa
ANNOTATION=genome.gtf

# [I commented out "Obtain the genome" section as the devlopment is just with the chromosome 22
# while the entire chromosome is too big. Uncomment the below if your machine
# has enough RAM (~100Gb) to work with the entire human genome] 

# Obtain the genome.
#wget -O .${GENOME_DIR}/${REFERENCE}.gz ${GENOME_LINK}
#gunzip -c .${GENOME_DIR}/${REFERENCE}.gz > .${GENOME_GENOME_DIRFOLDER}/${REFERENCE}
#rm .${GENOME_DIR}/${REFERENCE}.gz

# Obtain the genome annotation:
wget -O .${GENOME_DIR}/${ANNOTATION}.gz ${GENOME_ANNOTATION_LINK}
gunzip -c .${GENOME_DIR}/${ANNOTATION}.gz > .${GENOME_DIR}/${ANNOTATION}
rm .${GENOME_DIR}/${ANNOTATION}.gz

################

# 1) Create a simple fasta index file (.fai):
docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464 samtools \
    faidx ${REFERENCE}

# 2) Creaet a sequence dictionary used by GATK and Picard:
docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867 gatk CreateSequenceDictionary \
    -R ${REFERENCE} \
    -O ${REFERENCE%.*}.dict

# 3) Obtaining BED file from GTF/GFF (e.g., from Ensembl/NCBI).
# BED file describes genomic intervals (e.g., genes, exons, regions of interest)
docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/bedops:2.4.41--0451d22c61ea1547 bash -c "gff2bed \
< ${ANNOTATION} > ${ANNOTATION%.*}.bed"
