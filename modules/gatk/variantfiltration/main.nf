#!/usr/bin/env nextflow

params.outdir = "results/gatk/variantfiltration"

process GATK_VARIANT_FILTRATION {
    container "community.wave.seqera.io/library/bcftools_gatk4:d89f6490b3de65be"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(hc_variants), path(hc_variants_index)
        path genome
        path genome_fai
        path genome_dict

    output:
        tuple val(sample_id), path ("*.hc.pass.vcf.gz"), path ("*.hc.pass.vcf.gz.tbi"), emit: filtered_variants
        path ("*.hc.filt.vcf")

    script:
    """
    gatk VariantFiltration \
        -R ${genome} \
        -V ${hc_variants} \
        -O ${sample_id}.hc.filt.vcf \
        --filter-expression "QD < 2.0 || FS > 30.0 || MQ < 40.0" \
        --filter-name "RNAseqFilter" \
	    --verbosity ERROR

    bcftools filter -i 'FILTER="PASS"' -o ${sample_id}.hc.pass.vcf ${sample_id}.hc.filt.vcf

    bgzip ${sample_id}.hc.pass.vcf && tabix -p vcf ${sample_id}.hc.pass.vcf.gz
    """
}

