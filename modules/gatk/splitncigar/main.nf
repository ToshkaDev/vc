#!/usr/bin/env nextflow

params.outdir = "results/gatk/splitncigar"

/*
 * STAR
 */
 
process GATK_SPLIT_NCIGAR_READS {
    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(dupbam)
        path genome
        path genome_fai
        path genome_dict

    output:
        path "*", emit: cigar

    script:
    """
    gatk SplitNCigarReads \
	-R ${genome} \
	-I ${dupbam} \
	-O ${sample_id}.cigar.bam \
	--verbosity ERROR
    """
}

