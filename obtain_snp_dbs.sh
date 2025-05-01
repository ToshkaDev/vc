#!/bin/bash

SNP_DB_DIR=./data/snp_db
mkdir -p $SNP_DB_DIR

DB_SNP=Homo_sapiens_assembly38.dbsnp138.vcf.gz
STANDARD_INDELS=Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Following GATK Best Practices: Base Quality Score Recalibration (BQSR):
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890811--How-to-Recalibrate-base-quality-scores-run-BQSR
# known variants - dbSNP for known SNPs
DB_SNP_LINK=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/${DB_SNP}
# Mills and 1000G for known indels
STANDARD_INDELS_LINK=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/${STANDARD_INDELS}

# The following for the human genome will take ~3 min
# Obtain dbSNP v138 (1.45G)
wget -O ${SNP_DB_DIR}/${DB_SNP} ${DB_SNP_LINK}
# Obtain dbSNP v138 index file (2.21M)
wget -O ${SNP_DB_DIR}/${DB_SNP}.tbi ${DB_SNP_LINK}.tbi

# Obtain Mills and 1000G gold standard indels (19.73M)
wget -O ${SNP_DB_DIR}/${STANDARD_INDELS} ${STANDARD_INDELS_LINK}
# Mills and 1000G gold standard indels index file (1.43M)
wget -O ${SNP_DB_DIR}/${STANDARD_INDELS}.tbi ${STANDARD_INDELS_LINK}.tbi


