process ANNOTATE_VARIANTS {
    container 'bioconductor/bioconductor_docker:RELEASE_3_21'

    publishDir "results/annotatevariants", mode: 'copy'

    input:
        path file2sample
        path stats_folder
        path mut_filter_script
        path merged_stats

    output:
        path "*.tsv"

    script:
    """
    
    Rscript $mut_filter_script ${file2sample} ${stats_folder}
    """
}