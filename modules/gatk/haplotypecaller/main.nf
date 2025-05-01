#!/usr/bin/env nextflow

params.outdir = "results/gatk/haplotypecaller"

process GATK_HAPLOTYPE_CALLER {
    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(recalbam), path(recalbai)
        path genome
        path genome_fai
        path genome_dict
        path snpdb
        path snpdb_index

    output:
        path "*", emit: calledsnps

    script:
    """
	gatk HaplotypeCaller -R $genome \
    -I $recalbam \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp $snpdb \
    --dont-use-soft-clipped-bases \
    -O ${sample_id}.hc.vcf \
	--verbosity ERROR
    """
}

