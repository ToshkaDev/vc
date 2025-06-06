#!/bin/bash

# This script generates genome index according to the STAR manual.
# Please note that the genome index creation (especially of large genomes, such as the human)
# requires lots of RAM - This task together with the annotaions file was not be able to fininsh
# on a standard laptop with 30 GB RAM and 16 cores.

# Input file recommendations:
# FASTA: Ensembl "Primary Assembly" FASTA (e.g., Homo_sapiens.GRCh38.dna.primary_assembly.fa)
# GTF: Matching Ensembl GTF (e.g., Homo_sapiens.GRCh38.113.gtf)

# Set parameters GENOME, ANNOTATIONS (uncomment), and OVERHANG (uncomment) below 
# according to your genome, annotations file, and the length of your reads.
# --sjdbOverhang equals to the maximum read length you plan to map-1. 

if [ "$1" != "local" ] && [ "$1" != "hpc" ]; then
    echo "Specify either 'local' or 'hpc' as an argument to the script, ex., $0 local"
    exit 1
fi

DATA_DIR=/data
GENOME_DIR=/data/genome
GENOME=genome.fa
ANNOTATION=genome.gtf
OVERHANG=149

mkdir -p .${GENOME_DIR}/genome_index

echo "Creating genome STAR index ..."
# -w - setting working directory to ensure STAR writes logs and other temporary files 
# within the mounted volume
if [ "$1" == "local" ]; then
    echo "Docker found. Running with Docker..."
    docker run --rm -w ${GENOME_DIR}/genome_index -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/star:2.7.10b--90133b03b1960405 STAR \
        --runThreadN 8 \
        --runMode genomeGenerate \
        --genomeDir ${GENOME_DIR}/genome_index \
        --genomeFastaFiles ${GENOME_DIR}/${GENOME} \
        --sjdbGTFfile ${GENOME_DIR}/${ANNOTATION} \
        --sjdbOverhang ${OVERHANG}

elif [ "$1" == "hpc" ]; then
    echo "Singularity found. Running with Singularity..."
    singularity exec --pwd ${GENOME_DIR}/genome_index --bind .${DATA_DIR}:${DATA_DIR} docker://community.wave.seqera.io/library/star:2.7.10b--90133b03b1960405 STAR \
        --runThreadN 8 \
        --runMode genomeGenerate \
        --genomeDir ${GENOME_DIR}/genome_index \
        --genomeFastaFiles ${GENOME_DIR}/${GENOME} \
        --sjdbGTFfile ${GENOME_DIR}/${ANNOTATION} \
        --sjdbOverhang ${OVERHANG}
fi