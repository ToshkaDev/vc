process ANNOTATE_AND_ANALYZE_VARIANTS {
    container 'bioconductor/bioconductor_docker:RELEASE_3_21'

    publishDir "results/annotatevariants", mode: 'copy'

    input:
        path file2sample
        path stats_folder, dir: true

    output:
        path "tables/snv_indel_*.tsv"
        path "plots/*.pdf"
        path "tables/snv_indel_summary*.tsv"
        path "plots/hc_stats_library_size_vs_variant_number*.pdf"

    script:
    """
    Rscript scripts/analyze_variants.R ${file2sample} ${stats_folder}
    """
}