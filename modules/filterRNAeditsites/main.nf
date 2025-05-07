#!/usr/bin/env nextflow

params.outdir = "results/filterRNAeditsites"

process FILTER_RNA_EDIT_SITES {
    container "community.wave.seqera.io/library/tabix_vcftools:bfebada0dca23079"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(filtered_variants_gz), path(filtered_variants_gz_index)
        path rna_edit

    output:
        tuple val(sample_id), path ("*.hc.pass2.vcf.gz"), emit: rnaedit_filtered_vars

    script:
    """
    vcftools --gzvcf ${filtered_variants_gz} \
        --recode --exclude-positions $rna_edit --stdout | gzip -c > ${sample_id}.hc.pass2.vcf.gz
    """
}

