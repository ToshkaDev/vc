#!/usr/bin/env nextflow

params.outdir = "results/gatk/readgroups"

/*
 * STAR
 */
 
process GATK_ADD_REPLACE_READ_GROUPS {
    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path ("*.rg.bam"), path("*.rg.bai"), emit: rgs

    script:
    """
    gatk AddOrReplaceReadGroups \
    -I ${bam} \
    -O ${sample_id}.rg.bam \
    --SORT_ORDER coordinate \
    --RGID bc.eurofins.${sample_id} \
    --RGLB barstrand_specific \
    --RGPL illumina \
    --RGSM ${sample_id} \
    --RGPU flowcell.${sample_id} \
    --CREATE_INDEX True \
    --VERBOSITY ERROR
    """
}