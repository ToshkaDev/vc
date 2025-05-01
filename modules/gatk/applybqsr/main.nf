#!/usr/bin/env nextflow

params.outdir = "results/gatk/applybqsr"

/*
 * STAR
 */
 
process GATK_APPLY_BQSR {
    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(recaltable), path(cigarbam), path(cigarbai)
        path genome
        path genome_fai
        path genome_dict

    output:
        path "*", emit: appliedbqsr

    script:
    """
	gatk ApplyBQSR \
	--add-output-sam-program-record \
	-R $genome \
	-I $cigarbam  \
	--use-original-qualities \
	-O ${sample_id}.recal.bam \
	--bqsr-recal-file ${recaltable} \
	--verbosity ERROR
    """
}

