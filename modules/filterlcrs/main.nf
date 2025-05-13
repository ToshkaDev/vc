#!/usr/bin/env nextflow

params.outdir = "results/filterlcrs"

process FILTER_LCRS {
    container "community.wave.seqera.io/library/snpeff_snpsift:5.2--e2939680b6ff2466"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(rnaedit_filtered)
        path lcr_bed

    output:
        tuple val(sample_id), path("*.hc.pass2.lcr.vcf.gz"), emit: lcr_filtered_vars

    script:
    """    
    SnpSift intervals \
    -noLog -x -i ${rnaedit_filtered} \
    ${lcr_bed} | gzip > ${sample_id}.hc.pass2.lcr.vcf.gz
    """
}

