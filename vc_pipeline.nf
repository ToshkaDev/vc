#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input (file of input files, one per line)
params.input_csv = "${projectDir}/data/paired-end.csv"
params.genome_index="${projectDir}/data/genome/genome_index"

// Accessory files
// params.reference        = "${projectDir}/data/ref/ref.fasta"
// params.reference_index  = "${projectDir}/data/ref/ref.fasta.fai"
// params.reference_dict   = "${projectDir}/data/ref/ref.dict"
// params.intervals        = "${projectDir}/data/ref/intervals.bed"

include { FASTP } from './modules/fastp/main.nf'
include { STAR } from './modules/star/main.nf'


workflow {

    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [row.sample_id, file(row.fastq_1), file(row.fastq_2)] }

    FASTP(read_ch)

    STAR(FASTP.out.trimmed_reads, params.genome_index)
}
