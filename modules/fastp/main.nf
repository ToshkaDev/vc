#!/usr/bin/env nextflow

params.outdir = "results/fastp"

process FASTP {
    container "community.wave.seqera.io/library/fastp:0.24.1--6214360065b44e0b"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        path "*_fastqc.json", emit: json
        path "*_fastqc.html", emit: html
        tuple val(sample_id), path("*.R1.fq.gz"), path("*.R2.fq.gz"), emit: trimmed_reads

    script:
    """
    fastp -i ${read1} -I ${read2} \
        -o ${sample_id}.R1.fq.gz -O ${sample_id}.R2.fq.gz \
        -h ${sample_id}_fastqc.html -j ${sample_id}_fastqc.json \
        -q 20 -c -p --thread 4
    """
}