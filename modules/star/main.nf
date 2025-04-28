#!/usr/bin/env nextflow

/*
 * STAR
 */
process STAR {
    container "community.wave.seqera.io/library/star:2.7.10b--90133b03b1960405"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple path(read1), path(read2)

    output:
        path "*_fastqc.json", emit: json
        path "*_fastqc.html", emit: html
        tuple path("*.R1.fq.gz"), path("*.R2.fq.gz"), emit: trimmed_reads

    script:
    """
    fastp -i ${read1} -I ${read2} \
        -o ${read1.simpleName}.R1.fq.gz -O ${read1.simpleName}.R2.fq.gz \
        -h ${read1.simpleName}_fastqc.html -j ${read1.simpleName}_fastqc.json \
        -q 20 -c -p --thread 4
    """
}