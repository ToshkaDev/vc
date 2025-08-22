process ANNOTATE_VARIANTS {
    container 'bioliners/bioconductor_r3_21-vc:latest'

    publishDir "results/annotatevariants", mode: 'copy'

    input:
        path file2sample
        path stats_folder
        path mut_filter_script
        path merged_stats

    output:
        path "*.tsv", emit: annotated_vars

    script:
    """
    Rscript $mut_filter_script ${file2sample} ${stats_folder}
    """
}