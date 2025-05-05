#!/usr/bin/env nextflow

params.outdir = "results/star"

process STAR {
    container "community.wave.seqera.io/library/star:2.7.10b--90133b03b1960405"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)
        path genome_index


    output:
        path "*.tab"
        path "*Log.out"
        path "*Log.final.out", emit: logs
        tuple val(sample_id), path("*.bam"), emit: bam

    script:
    """
     STAR \
    --genomeDir ${genome_index} \
    --runThreadN 14 \
    --readFilesIn ${read1} ${read2} \
    --readFilesCommand zcat \
    --sjdbOverhang 149 \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 20000000000 \
    --outSAMunmapped Within \
    --outFileNamePrefix "${sample_id}_"

    """
}