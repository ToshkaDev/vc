process COLLECT_STATS {

    publishDir "results/collectstats", mode: 'copy'

    input:
        tuple val(sample_id), path ("*.hc.vcf"), path(_)
        
        tuple val(sample_id), path(hc_pass_vcf), path(_), path(hc_pass2_vcf), path(hc_pass2_lcr_vcf), 
            path(hc_pass2_lcr_1kG_vcf), path(hc_pass2_lcr_dbsnp_vcf)

        tuple val(sample_id), path(pass_vcf_maf), path(pass2_vcf_maf), path(pass2_lcr_vcf_maf), 
            path(pass2_lcr_1kG_vcf_maf), path(pass2_lcr_dbsnp_vcf_maf)

        tuple val(sample_id), path(log_files)

    output:
        path("*")

    script:
    """
    # For log files
    for f in ${log_files[*]}; do
        echo "Processing: \$f"
        # do something with \$f
    done

    """
}