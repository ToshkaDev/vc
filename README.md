# 🧬 Nextflow Variant Calling Pipeline

# Status: 🚧 Heavy development underway — expect frequent updates and improvements.

This repository contains a modular, reproducible Nextflow pipeline for variant discovery from high-throughput sequencing data. The pipeline follows the Broad Institute’s GATK Best Practices but is designed to be flexible, allowing use of alternative tools like DeepVariant for variant calling. It processes raw FASTQ files through to variant calling and includes extensive quality control, alignment, base recalibration, and variant calling steps — with testing integrated at every stage.

## 📋 Features

Reproducible and scalable

GATK Best Practices-compatible, but supports alternative variant callers (e.g., DeepVariant)

Includes QC, alignment, duplicate marking, BQSR, and variant calling

Per-step testing with small datasets

Compatible with Docker, Singularity, and Conda

Ready for use on local, HPC, or cloud platforms

## System Requirements

## Hardware requirements

- CPU: 8+ cores recommended (parallel execution supported via Nextflow)
- Memory:
  - Minimum: 32 GB RAM
  - Recommended:
    - 64–128 GB RAM for standard variant calling with DeepVariant or GATK
    - >150 GB RAM required for STAR genome indexing with annotation file and large-scale RNA-seq alignment (e.g., human hg38)
- Disk: Minimum 100 GB free disk space per WGS sample
- GPU (optional): Recommended for DeepVariant acceleration (NVIDIA GPU with CUDA support)

> ⚠️ Requirements depend on dataset size and selected tools. Cloud/HPC deployment is recommended for high-throughput analyses.


## Software requirements

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
```./create_required_files.sh``` obtains the human genome and its annotation, creates fai index, and sequence dictionary used by GATK and Picard, and crates a BED file from the provided annotation file. The script can work with other genomes if corresponding links are provided: options GENOME_LINK and GENOME_ANNOTATION_LINK.

```./create_genome_index.sh``` creates a STAR genome index, with chr and scaffolds, primary assembly, and transcriptome as recommended by the STAR manual.

```./obtain_snp_dbs.sh``` obtains known variants (dbSNP for known SNPs) and known indels (Mills and 1000G) following [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811--How-to-Recalibrate-base-quality-scores-run-BQSR
)

The data folder currently includes human chr22. To work with the entire human genome uncomment the 'Obtain the genome' before running ```./create_genome_index.sh```.

**Now run the pipeline**:
```
nextflow run vc_pipeline.nf
```

Supports execution on:

Local workstation

HPC (SLURM/LSF) (coming soon)

Cloud (AWS, Google Cloud) (coming soon)

## 🛠 Workflow Overview

* Quality Control & Trimming (fastp)

* Aggregated QC Reporting (MultiQC)

* Alignment (BWA for DNA, STAR for RNA)

* Post-processing (Sort)

* Pre-processing for variatn calling (AddOrReplaceReadGroups, MarkDuplicates, SplitNCigarReads (RNA only))

* Base Quality Score Recalibration (BaseRecalibrator, ApplyBQSR) (optional if using DeepVariant)

* Variant Calling

Options: GATK HaplotypeCaller, DeepVariant (coming soon)

## 📦 Inputs

FASTQ files (paired-end; single-end - coming soon)

Reference genome (e.g. Homo_sapiens_assembly38.fasta)

Known variant sites: dbSNP, Mills_and_1000G_gold_standard.indels

Optional sample sheet CSV with sample metadata

## 📤 Outputs

Quality metrics

Aligned, deduplicated, and recalibrated BAM files

Raw and recalibrated VCF/GVCF variant calls

Log files and structured QC reports

## ✅ Testing

Each module has small, testable datasets

Continuous validation via nf-test

Easy test runs for new environments or updates

## 🧱 Built With

Nextflow

shell

GATK 4.5+

DeepVariant

Docker

## 📚 Documentation

docs/ folder (coming soon)

Example input files and configuration templates

## 🙌 Acknowledgements

Broad Institute GATK team

Google DeepVariant team

nf-core community inspiration

## 📌 License

This pipeline is open-source and available under the MIT License.

For questions, suggestions, or to contribute: please open an issue or submit a pull request.
