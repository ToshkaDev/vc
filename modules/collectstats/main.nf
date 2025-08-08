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
        tuple val(sample_id), path ("*.hc.vcf")
        tuple val(sample_id), path(hc_pass_vcf), path(_), path(hc_pass2_vcf), path(hc_pass2_lcr_vcf), 
            path(hc_pass2_lcr_1kG_vcf), path(hc_pass2_lcr_dbsnp_vcf)
        tuple val(sample_id), path(pass_vcf_maf), path(pass2_vcf_maf), path(pass2_lcr_vcf_maf), 
            path(pass2_lcr_1kG_vcf_maf), path(pass2_lcr_dbsnp_vcf_maf)
        path("*.txt")

    script:
    """

    # For log files
    for f in ${log_files[*]}; do
        echo "Processing: \$f"
        # do something with \$f
    done

    # q20
    zgrep -c "^chr" *.hc.vcf.gz > stats_hc_all.txt
    # q20 + depth
    zgrep -c "^chr" *.hc.pass.vcf.gz > stats_hc_pass.txt
    # q20 + depth + rna_edit
    zgrep -c "^chr" *.hc.pass2.vcf.gz > stats_hc_pass2.txt
    # q20 + depth + rna_edit + lcr
    zgrep -c "^chr" *.hc.pass2.lcr.vcf.gz > stats_hc_lcr.txt
    # q20 + depth + rna_edit + lcr + 1000G
    zgrep -c "^chr" *.hc.pass2.lcr.1kG.vcf.gz > stats_hc_1kG.txt
    # q20 + depth + rna_edit + lcr + dbsnp
    zgrep -c "^chr" *.hc.pass2.lcr.dbsnp.vcf.gz > stats_hc_dbsnp.txt

    """
}