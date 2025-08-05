#!/bin/bash

if [ "$1" != "local" ] && [ "$1" != "hpc" ]; then
    echo "Specify either 'local' or 'hpc' as an argument to the script, ex., $0 local"
    exit 1
fi

# Broad Institute's compatable files:
GENOME_LINK=https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
GENOME_ANNOTATION_LINK=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz

DATA_DIR=/data
GENOME_DIR=/data/genome
REFERENCE=genome.fa
ANNOTATION=genome.gtf

# The devlopment is just with the chromosome 22. Empty the .${GENOME_DIR}/${REFERENCE} folder before running the script if you want to work with the entire human genome 

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
echo "Creating simple fasta index, sequence, and BED files ..."

if [ "$1" == "local" ] && command -v docker &> /dev/null; then
    echo "Docker found. Running with Docker..."
    # 1) Create a simple fasta index file (.fai):
    echo "Creating simple fasta index file ..."
    docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464 samtools \
        faidx ${REFERENCE}

    # 2) Creaet a sequence dictionary used by GATK:
    echo "Creating sequence dictionary ..."
    docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867 gatk CreateSequenceDictionary \
        -R ${REFERENCE} \
        -O ${REFERENCE%.*}.dict

    # 3) Obtaining BED file from GTF/GFF (e.g., from Ensembl/NCBI).
    # BED file describes genomic intervals (e.g., genes, exons, regions of interest)
    echo "Creating BED file (describes genomic intervals) ..."
    docker run --rm -w ${GENOME_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/bedops:2.4.41--0451d22c61ea1547 bash -c "gff2bed \
    < ${ANNOTATION} > ${ANNOTATION%.*}.bed"

elif [ "$1" == "hpc" ] && command -v singularity &> /dev/null; then
    echo "Singularity found. Running with Singularity..."
    echo "Creating simple fasta index file ..."
    singularity exec --pwd ${GENOME_DIR} --bind .${DATA_DIR}:${DATA_DIR} docker://community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464 samtools \
        faidx ${REFERENCE}

    echo "Creating sequence dictionary ..."
    singularity exec --pwd ${GENOME_DIR} --bind .${DATA_DIR}:${DATA_DIR} docker://community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867 gatk CreateSequenceDictionary \
        -R ${REFERENCE} \
        -O ${REFERENCE%.*}.dict

    echo "Creating BED file (describes genomic intervals) ..."
    singularity exec --pwd ${GENOME_DIR} --bind .${DATA_DIR}:${DATA_DIR} docker://community.wave.seqera.io/library/bedops:2.4.41--0451d22c61ea1547 bash -c "gff2bed \
    < ${ANNOTATION} > ${ANNOTATION%.*}.bed"
else
    echo "No supported containerization tool matching argument '$1' was found in the environment."
fi

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

################
# Download VEP cache. The VEP cache contains precompiled gene annotation data that 
# VEP (Variant Effect Predictor) uses to rapidly annotate variants offline, without repeatedly querying remote Ensembl databases.
# The VEP cache includes transcript annotations and gene models.
VEP_CACHE=.${DATA_DIR}/vep_cache

mkdir -p $VEP_CACHE
echo "Downloading and staging VEP cache with precompiled gene annotation and gene models. Depending on the internet connection this
may take up to several hours (~14 GB volume of compressed data)"
if [ "$1" == "local" ] && command -v docker &> /dev/null; then
    echo "Docker found. Running with Docker ..."
    # vep_install script inside the ensembl-vep_vcf2maf image is used
    docker run --rm -v $VEP_CACHE:/vep_cache community.wave.seqera.io/library/ensembl-vep_vcf2maf:6836d39d9a125f9f bash \
        -c "export VEP_DOWNLOAD_METHOD=http && yes n | /opt/conda/bin/vep_install -a cf -s homo_sapiens -y GRCh38 -c /vep_cache"
elif [ "$1" == "hpc" ] && command -v singularity &> /dev/null; then
    echo "Singularity found. Running with Singularity ..."
    singularity exec --bind $VEP_CACHE:/vep_cache docker://community.wave.seqera.io/library/ensembl-vep_vcf2maf:6836d39d9a125f9f bash \
        -c "export VEP_DOWNLOAD_METHOD=http && yes n | /opt/conda/bin/vep_install -a cf -s homo_sapiens -y GRCh38 -c /vep_cache"
else
    echo "No supported containerization tool matching argument '$1' was found in the environment."
fi