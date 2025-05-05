#!/bin/bash

SNP_DB_DIR=./data/snp_db
mkdir -p $SNP_DB_DIR

DB_SNP_NAME="snpdb.vcf.gz"
INDELS_NAMES="indels.vcf.gz"

# Following GATK Best Practices: Base Quality Score Recalibration (BQSR):
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890811--How-to-Recalibrate-base-quality-scores-run-BQSR
# known variants - dbSNP for known SNPs
DB_SNP_LINK=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
# Mills and 1000G for known indels
STANDARD_INDELS_LINK=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# The following for the human genome will take ~3 min
# Obtain dbSNP v138 (1.45G)
echo "Downloading dbSNP v138 and its index file (1.45G + 2.21M) ..."
wget -O ${SNP_DB_DIR}/${DB_SNP_NAME} ${DB_SNP_LINK}
# Obtain dbSNP v138 index file (2.21M)
wget -O ${SNP_DB_DIR}/${DB_SNP_NAME}.tbi ${DB_SNP_LINK}.tbi

# Obtain Mills and 1000G gold standard indels (19.73M)
echo "Downloading Mills and 1000G gold standard indels and its index file (9.73M + 1.43M) ..."
wget -O ${SNP_DB_DIR}/${INDELS_NAMES} ${STANDARD_INDELS_LINK}
# Mills and 1000G gold standard indels index file (1.43M)
wget -O ${SNP_DB_DIR}/${INDELS_NAMES}.tbi ${STANDARD_INDELS_LINK}.tbi


