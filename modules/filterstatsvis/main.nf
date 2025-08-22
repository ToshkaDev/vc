process FILTER_STATS_VISUALIZE {
    container 'bioliners/bioconductor_r3_21-vc:latest'

    publishDir "results/filterstatsvis", mode: 'copy'

    input:
        path stats_folder
        path annotatevars_folder
        path starstats_folder
        path filter_stats_vis_script
        path annotated_vars


    output:
        path "*.pdf", optional: true

    script:
    """
    Rscript $filter_stats_vis_script ${stats_folder} ${annotatevars_folder} ${starstats_folder}

    """
}