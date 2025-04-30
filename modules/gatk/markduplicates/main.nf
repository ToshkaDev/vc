#!/usr/bin/env nextflow

params.outdir = "results/gatk/markduplicates"

/*
 * STAR
 */
 
process GATK_MARK_DUPLICATES {
    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(bamrgs), path(bairgs)

    output:
        path("*_metrics.txt")
        tuple val(sample_id), path("*_dup.bam"), emit: dups

    script:
    """
    gatk MarkDuplicates \
    -I ${bamrgs} \
    -O ${sample_id}.mark_dup.bam \
    -M ${sample_id}.mark_dup_metrics.txt \
    --VERBOSITY ERROR
    """
}

