#!/usr/bin/env nextflow

params.outdir = "results/gatk/baserecalibrate"

/*
 * STAR
 */
 
process GATK_BASE_RECALIBRATOR {
    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(cigarbam), path(cigarbai)
        path genome
        path genome_fai
        path genome_dict
        path snpdb
        path snpdb_index
        path indel
        path indel_index

    output:
        tuple val(sample_id), path ("*.recal.table"), path(cigarbam), path(cigarbai), emit: baserecal

    script:
    """
	gatk BaseRecalibrator \
	-R $genome \
	-I $cigarbam \
	--use-original-qualities \
	--known-sites $snpdb \
	--known-sites $indel \
	-O ${sample_id}.recal.table \
	--verbosity ERROR
    """
}

