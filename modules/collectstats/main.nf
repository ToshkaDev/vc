process ANNOTATE_AND_ANALYZE_VARIANTS {

    publishDir "results/collectstats", mode: 'copy'

    input:
        tuple val(sample_id), path(hc_pass_vcf), path(_), path(hc_pass2_vcf), path(hc_pass2_lcr_vcf), 
        path(hc_pass2_lcr_1kG_vcf), path(hc_pass2_lcr_dbsnp_vcf) 
    
    output:


    script:
    """
    """
}