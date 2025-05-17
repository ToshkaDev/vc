#!/usr/bin/env nextflow

params.outdir = "results/filtercommonvars"

process FILTER_COMMON_VARIANTS {
    container "community.wave.seqera.io/library/tabix_vcftools:bfebada0dca23079"

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(lcr_filtered)
        path(gpp3snp_af01_filter)
        path(ncbi_dbsnp_filter)

    output:
        tuple val(sample_id), path("*.1kG.snp.vcf.gz"), path("*.1kG.vcf.gz"), path("*.dbsnp.snp.vcf.gz"),  path("*.dbsnp.vcf.gz"), emit: common_filtered_vars

    script:
    """
    # filter out common snps from called snvs, without lcr, without rna-edit sites:

    # without indels for mutsigs, gatk 1000 genomes
    vcftools --gzvcf $lcr_filtered --exclude-positions ${gpp3snp_af01_filter} --remove-indels \
        --recode --recode-INFO-all --stdout | bgzip > ${sample_id}.hc.pass2.lcr.1kG.snp.vcf.gz
    
    # with indels for waterfall, gatk 1000 genomes
    vcftools --gzvcf $lcr_filtered --exclude-positions ${gpp3snp_af01_filter}\
        --recode --recode-INFO-all --stdout | bgzip > ${sample_id}.hc.pass2.lcr.1kG.vcf.gz

    # without indels for waterfall, nih dbsnp
    vcftools --gzvcf $lcr_filtered --exclude-positions ${ncbi_dbsnp_filter} --remove-indels \
        --recode --recode-INFO-all --stdout | bgzip > ${sample_id}.hc.pass2.lcr.dbsnp.snp.vcf.gz

    # with indels for waterfall, nih dbsnp
    vcftools --gzvcf $lcr_filtered --exclude-positions ${ncbi_dbsnp_filter} \
        --recode --recode-INFO-all --stdout | bgzip  > ${sample_id}.hc.pass2.lcr.dbsnp.vcf.gz
    """
}

