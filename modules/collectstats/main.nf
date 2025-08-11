process COLLECT_STATS {

    publishDir "results/collectstats", mode: 'copy'

    input:
        path(star_stats)
        tuple val(sample_id), path(hc_vcf), path(_)
        tuple val(sample_id), path(hc_pass_vcf), path(_), path(hc_pass2_vcf), path(hc_pass2_lcr_vcf), 
            path(hc_pass2_lcr_1kG_vcf), path(hc_pass2_lcr_dbsnp_vcf)
        tuple val(sample_id), path(pass_vcf_maf), path(pass2_vcf_maf), path(pass2_lcr_vcf_maf), 
            path(pass2_lcr_1kG_vcf_maf), path(pass2_lcr_dbsnp_vcf_maf)
        tuple val(sample_id), path(log_files)

    output:
        path(star_stats)
        tuple val(sample_id), path(hc_vcf)
        tuple val(sample_id), path(hc_pass_vcf), path(_), path(hc_pass2_vcf), path(hc_pass2_lcr_vcf), 
            path(hc_pass2_lcr_1kG_vcf), path(hc_pass2_lcr_dbsnp_vcf)
        tuple val(sample_id), path(pass_vcf_maf), path(pass2_vcf_maf), path(pass2_lcr_vcf_maf), 
            path(pass2_lcr_1kG_vcf_maf), path(pass2_lcr_dbsnp_vcf_maf)
        path "${sample_id}_*.txt", emit: stats_files

    script:
    """
    set -euo pipefail
    # q20
    printf "%s:%s\\n" "$hc_vcf" "\$(zgrep -c '^chr' ${hc_vcf} || echo 0)" > ${sample_id}_hc_all.txt
    # q20 + depth
    printf "%s:%s\\n" "$hc_pass_vcf" "\$(zgrep -c '^chr' ${hc_pass_vcf} || echo 0)" > ${sample_id}_hc_pass.txt
    # q20 + depth + rna_edit
    printf "%s:%s\\n" "$hc_pass2_vcf" "\$(zgrep -c '^chr' ${hc_pass2_vcf} || echo 0)" > ${sample_id}_hc_pass2.txt
    # q20 + depth + rna_edit + lcr
    printf "%s:%s\\n" "$hc_pass2_lcr_vcf" "\$(zgrep -c '^chr' ${hc_pass2_lcr_vcf} || echo 0)" > ${sample_id}_hc_lcr.txt
    # q20 + depth + rna_edit + lcr + 1000G
    printf "%s:%s\\n" "$hc_pass2_lcr_1kG_vcf" "\$(zgrep -c '^chr' ${hc_pass2_lcr_1kG_vcf} || echo 0)" > ${sample_id}_lcr_1kG.txt
    # q20 + depth + rna_edit + lcr + dbsnp
    printf "%s:%s\\n" "$hc_pass2_lcr_dbsnp_vcf" "\$(zgrep -c '^chr' ${hc_pass2_lcr_dbsnp_vcf} || echo 0)" > ${sample_id}_lcr_dbsnp.txt
    """
}

process MERGE_STATS {
    publishDir "results/collectstats", mode: 'copy'
        
    input:
        path(stats)

    output:
        path "stats_*.txt"

    script:
    """
    # q20
    cat *_hc_all.txt  > stats_hc_all.txt
    # q20 + depth
    cat *_hc_pass.txt > stats_hc_pass.txt
    # q20 + depth + rna_edit
    cat *_hc_pass2.txt > stats_hc_pass2.txt
    # q20 + depth + rna_edit + lcr
    cat *_hc_lcr.txt > stats_hc_lcr.txt
    # q20 + depth + rna_edit + lcr + 1000G
    cat *_lcr_1kG.txt > stats_hc_1kG.txt
    # q20 + depth + rna_edit + lcr + dbsnp
    cat *_lcr_dbsnp.txt > stats_hc_dbsnp.txt
    """

}