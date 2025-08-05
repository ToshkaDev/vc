process ANNOTATE_AND_ANALYZE_VARIANTS {
    container 'bioconductor/bioconductor_docker:RELEASE_3_21'

    publishDir "results/annotatevariants", mode: 'copy'

    input:
        path "snp/*.vcf.gz"
        path "snp/*.vcf.maf.gz"
        path "file2sample.csv"
        path "tables/star_stat_summary.tsv"
        path "snp/stats_hc_*.txt"
    
    output:
        path "tables/snv_indel_*.tsv"
        path "plots/*.pdf"
        path "tables/snv_indel_summary*.tsv"
        path "plots/hc_stats_library_size_vs_variant_number*.pdf"

    script:
    """
    Rscript scripts/analyze_variants.R
    """
}