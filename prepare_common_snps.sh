#!/bin/bash

# This script Prepares a file containing the mapping of RefSeq genome assembly identifiers to chromosome names and downloads and prepares 
# the 1000 Genomes Project Phase 3 data variant sites and the latest NCBI dbsnp set (~30 GB; download time ~45 min).
# These files will be used to filter out common snps

if [ "$1" != "docker" ] && [ "$1" != "singularity" ]; then
    echo "Specify either 'singularity' or 'docker' as an argument to the script, ex., $0 docker"
    exit 1
fi

DATA_DIR=/data
# Download the assembly report file to prepare a file containing the mapping of RefSeq genome assembly identifiers (e.g. GCF attachments) to chromosome names
# To be used for filtering common snps from called variants
ASSEMBLY_DIR=.${DATA_DIR}/assembly_mapping
ASSEMBLY_VERSION=GCF_000001405.40_GRCh38.p14
ASSEMBLY_REPORT_LINK=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/${ASSEMBLY_VERSION}/${ASSEMBLY_VERSION}_assembly_report.txt
ASSEMBLY_REPORT_FILE=assembly_report_file.tsv
MAPPING_FILE=gcf_refseq2chr.tsv

mkdir -p ${ASSEMBLY_DIR}

if [ ! -f ${ASSEMBLY_DIR}/${MAPPING_FILE} ]; then
    echo "Download the assembly report file to prepare a file containing the mapping of RefSeq genome assembly ids to chromosome names ..."
    wget -O ${ASSEMBLY_DIR}/${ASSEMBLY_REPORT_FILE} ${ASSEMBLY_REPORT_LINK}
    awk -F '\t' '$0 !~ /^#/ && $5 != "na" {print $7 "\t" $1}' ${ASSEMBLY_DIR}/${ASSEMBLY_REPORT_FILE} > ${ASSEMBLY_DIR}/${MAPPING_FILE}
    rm ${ASSEMBLY_DIR}/${ASSEMBLY_REPORT_FILE}
fi

################
# Downloading the latest NCBI dbsnp file for filtering common snps from called variants
DBSNP_DIR=.${DATA_DIR}/ncbi_dbsnp
DBSNP_LINK=https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
DBSNP_FILE=ncbi_dbsnp_latest.gz
FINAL_FILE=ncbi_dbsnp_latest.common.chr.tsv

mkdir -p ${DBSNP_DIR}

# If the final file is not available perform the following steps
if [ ! -f ${DBSNP_DIR}/${FINAL_FILE} ]; then
    # First check if the compressed prepared files are available and uncompress and concatenate them to the final file if they are
    if [ -f ${DBSNP_DIR}/${FINAL_FILE}.1.gz ] && [ -f ${DBSNP_DIR}/${FINAL_FILE}.2.gz ]; then
        echo "The ${DBSNP_DIR}/${FINAL_FILE}.1.gz and ${DBSNP_DIR}/${FINAL_FILE}.2.gz file are detected. Unpacking and concatenating ..."
        gunzip -c ${DBSNP_DIR}/${FINAL_FILE}.1.gz > ${DBSNP_DIR}/${FINAL_FILE}
        gunzip -c ${DBSNP_DIR}/${FINAL_FILE}.2.gz >> ${DBSNP_DIR}/${FINAL_FILE}
    # If those files are not available then perform all the necessary steps to obtain and prepare the final file
    else
        echo "Downloading gz compressed NBCI dbsnp database file (> 27.5 G, 10-20 min) ..."
        wget -O ${DBSNP_DIR}/${DBSNP_FILE} ${DBSNP_LINK}

        echo "Prepare file for filtering called variants (~10-15 min; requires RAM; can be helped through periodic cache dropping) ..."
        zgrep "^#" ${DBSNP_DIR}/${DBSNP_FILE} > ${DBSNP_DIR}/${DBSNP_FILE%.*}.common.vcf

        echo "Label COMMON for variants with AF>=0.01 (1%) (10-15 min; requires RAM; can be helped through periodic cache dropping) ... "
        zgrep ";COMMON" ${DBSNP_DIR}/${DBSNP_FILE} >> ${DBSNP_DIR}/${DBSNP_FILE%.*}.common.vcf

        # get chr and position for common variants (>0.01)
        refseq=$(cut -f1 ${ASSEMBLY_DIR}/${MAPPING_FILE})
        IFS=' ' read -ra refseq <<< $(echo $refseq)
        chr=$(cut -f2 ${ASSEMBLY_DIR}/${MAPPING_FILE})
        IFS=$' ' read -ra chr <<< $(echo $chr)

        touch ${DBSNP_DIR}/${DBSNP_FILE%.*}.common.chr.tsv

        echo "Prpare a file with common variants replacing refseq ids with corresponding chromosome ids in the NCBI dbsnp file (~40 min) ..."
        for i in ${!refseq[@]}; do
            date
            echo ${chr[$i]}
            grep ${refseq[$i]} ${DBSNP_DIR}/${DBSNP_FILE%.*}.common.vcf | sed "s/${refseq[$i]}/${chr[$i]}/g" | cut -f1,2 >> ${DBSNP_DIR}/${DBSNP_FILE%.*}.common.chr.tsv
        done
    fi
fi

#####
# Downloading a variant call format (VCF) file derived from the 1000 Genomes Project Phase 3 data, 
# aligned to the GRCh38 (hg38) human genome reference. This file includes variant sites with an allele frequency (AF) 
# of at least 1% (AF ≥ 0.01). 
# The final prepared file with variants with allele frequency (AF) ≥ 0.01 will be used for variant filtering 
GPP_DIR=${DATA_DIR}/gpp_vcf
GPP_LINK=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase3_v4_20130502.sites.hg38.vcf
GPP_FILE=1000G_phase3_sites.hg38.vcf

mkdir -p .${GPP_DIR}

if [ ! -f .${GPP_DIR}/${GPP_FILE} ]; then
    echo "Downloading the VCF file derived from the 1000 Genomes Project Phase 3 data and its index (5.3G + 1.09M; ~4-5 min) ..."
    wget -O .${GPP_DIR}/${GPP_FILE} ${GPP_LINK}
    wget -O .${GPP_DIR}/${GPP_FILE}.idx ${GPP_LINK}.idx

    echo "Filter variants with allele frequency (AF) ≥ 0.01 to be used for variant filtering ..."
    if [ "$1" == "docker" ]; then
        echo "Docker found. Running with Docker..."
        docker run --rm -w ${GPP_DIR} -v .${DATA_DIR}:${DATA_DIR} community.wave.seqera.io/library/bcftools_gatk4:d89f6490b3de65be \
            bcftools view -i 'INFO/AF>=0.01' ${GPP_DIR}/${GPP_FILE} -o ${GPP_DIR}/${GPP_FILE%.*}.af01.vcf

    elif [ "$1" == "singularity" ]; then
        echo "Singularity found. Running with Singularity..."
        singularity exec --pwd ${GPP_DIR} --bind .${DATA_DIR}:${DATA_DIR} docker://community.wave.seqera.io/library/bcftools_gatk4:d89f6490b3de65be \
            bcftools view -i 'INFO/AF>=0.01' ${GPP_DIR}/${GPP_FILE} -o ${GPP_DIR}/${GPP_FILE%.*}.af01.vcf
    fi
fi