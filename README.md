# 🧬 Nextflow Variant Calling Pipeline
<!-- Uncomment when GitHub Actions CI is configured
![Build Status](https://img.shields.io/github/actions/workflow/status/ToshkaDev/vc/vc_pipeline.nf?branch=main)
-->
![Nextflow](https://img.shields.io/badge/Nextflow-%E2%9C%94%20v24.10%2B-brightgreen)
![License](https://img.shields.io/github/license/ToshkaDev/vc)
![System Requirements](https://img.shields.io/badge/system-Java%2011%2B%20%7C%20Linux%2FmacOS%20%7C%20Nextflow%2024.10%2B-blue)
![Status](https://img.shields.io/badge/status-active%20development-yellow)
![GATK Support](https://img.shields.io/badge/GATK-✓%20Supported-blueviolet)
<!-- Uncomment when ready
![Docker](https://img.shields.io/badge/container-Docker%20%7C%20Singularity-orange)
-->
<!-- Uncomment when ready
![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)
-->
<!-- Uncomment when implemented
![DeepVariant Support](https://img.shields.io/badge/DeepVariant-✓%20Supported-green)
-->

# Status: Heavy development underway — expect frequent updates and improvements.

This repository contains a modular, reproducible Nextflow pipeline for variant discovery and analysis from high-throughput sequencing data. The pipeline follows the Broad Institute’s GATK Best Practices. It processes raw FASTQ files through variant calling and includes extensive quality control, alignment, preprocessing, base recalibration, variant calling, and subsequent filtering steps — with testing integrated at every stage.

## 📋 Features

Reproducible and scalable

GATK Best Practices-compatible

Includes QC, alignment, duplicate marking, BQSR, variant calling, and multi-stage result filtering

Per-step testing with small datasets

Ready for use on local, HPC, or cloud platforms

## System Requirements

### Hardware requirements

- CPU: 8+ cores recommended (parallel execution supported via Nextflow)
- Memory:
  - Minimum: 32 GB RAM
  - Recommended:
    - 64–128 GB RAM for standard variant calling withGATK
    - >150 GB RAM required for STAR genome indexing with annotation file and large-scale RNA-seq alignment (e.g., human hg38)
- Disk: Minimum 100 GB free disk space
<!--- GPU (optional): Recommended for DeepVariant acceleration (NVIDIA GPU with CUDA support) -->

> ⚠️ Requirements depend on dataset size and selected tools. Cloud/HPC deployment is recommended for high-throughput analyses.


### Software requirements

Nextflow 24.10+

Java 11+

Docker
 
## 🚀 Quick Start
```
git clone https://github.com/ToshkaDev/vc.git

cd vc

./create_required_files.sh

./create_genome_index.sh

./obtain_snp_dbs.sh
```

`./create_required_files.sh`:
- obtains the human genome and its annotation
- creates fai index and sequence dictionary used by GATK (and Picard)
- crates a BED file from the provided annotation file. The script can work with other genomes if corresponding links are provided (variables GENOME_LINK and GENOME_ANNOTATION_LINK)
- downloads and prepares an RNA edit sites file for the human genome. For other genomes provide an appropriate link (the RNA_EDIT_SITES_LINK variable)
- downloads and prepares a low complexity regions file for the human genome. For other genomes provide and appropriate link (the LCR_LINK variable)
- prepares a file containing the mapping of RefSeq genome assembly identifiers to chromosome names (will be used to prepare the NCBI dbsnp set for filtering common snps)
- downloads and prepares the 1000 Genomes Project Phase 3 data variant sites and the latest NCBI dbsnp set (both will be used to filter out common snps)

`./create_genome_index.sh` creates a STAR genome index, with chr and scaffolds, primary assembly, and transcriptome as recommended by the STAR manual.

The **data/genome** folder currently includes human chr22. To work with the entire human genome simply delete the genome.fa file in the ./data/genome folder before running `create_required_files.sh`. Genome indexing with STAR for the entire human geome (GRCh38 + GENCODE GTF annotaion file) will require ~100-150 GB of RAM. On a machine with 32 threads, SSD, --sjdbOverhang 100 (common for 101 bp reads), and 128–150 GB RAM available this can take ~1.5–2 hours.

`./obtain_snp_dbs.sh` 
- obtains known variants (dbSNP for known SNPs)
- known indels (Mills and 1000G) following [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811--How-to-Recalibrate-base-quality-scores-run-BQSR
)

**Once the preparatory steps are complete, start the pipeline**:
```
nextflow run vc_pipeline.nf
```

## 🛠 Workflow Overview

- Quality Control & Trimming (fastp)

- Aggregated QC Reporting (MultiQC)

- Alignment (STAR for RNA)

- Pre-processing for variant calling:
  - AddOrReplaceReadGroups
  - MarkDuplicates
  - SplitNCigarReads
  - Base Quality Score Recalibration (BaseRecalibrator, ApplyBQSR)
  
- Variant Calling (GATK HaplotypeCaller)
- Called vairant filtration:
  - Filtering low-quality variants (GATK VariantFiltration)
  - Filtering RNA edit sites (using vcftools and a set of RNA edit sites)
  - Filtering low complexity regions (using SnpSift and a corresponding set of marked low complexity regions)
  - Filtering common snps (using the 1000 Genomes Project Phase 3 data variant sites and the latest NCBI dbsnp set)

## 📦 Inputs

The input data listed below is obtained/prepared for the human genome by simply running the setup scripts provided in this repository. The scripts can extract and prepare all the necessary files for other genomes if the appropriate links are provided in these scripts.

FASTQ files (paired-end; single-end - coming soon)

Reference genome (e.g. Homo_sapiens_assembly38.fasta)

Known variant sites: dbSNP, Mills_and_1000G_gold_standard.indels

RNA edits sites (basic knowledge about RNA-edit sites can be found here: https://academic.oup.com/nar/article/49/D1/D1012/5940507?login=true)

Low complexity regions (information about LCR can be found here: https://academic.oup.com/bioinformatics/article/30/20/2843/2422145?login=true)

Sample sheet CSV with sample metadata

## 📤 Outputs

Quality metrics

Aligned, deduplicated, and recalibrated BAM files

Raw and recalibrated VCF/GVCF variant calls

Log files and structured QC reports

## ✅ Testing

Each module is equipped with a test section to ensure continuous validation. To be able to run tests install [nf-test](https://www.nf-test.com/installation/) first

Per module test example: ``` nf-test test modules/gatk/readgroups/tests/main.nf.test ```

The workflow-level testing: ``` nf-test test tests/vc_pipeline.nf.test ```

Running all the tests at once: ``` nf-test test ```

## Built With

Nextflow, shell, Docker; GATK 4.5+, fastp, SnpSift, vcftools, STAR, MultiQC

## Documentation

docs/ folder (coming soon)

Example input files and configuration templates

## Acknowledgements

Broad Institute GATK team

nf-core community inspiration

Developers of various beautiful tools used in variant calling and analysis

## License

This pipeline is open-source and available under the MIT License.

For questions, suggestions, or to contribute: please open an issue or submit a pull request.
