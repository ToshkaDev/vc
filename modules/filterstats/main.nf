process FILTER_STATS {
    container 'bioliners/bioconductor_r3_21-vc:latest'

    publishDir "results/filterstats", mode: 'copy'

    input:
        path stats_folder
        path annotatevars_folder
        path filter_stats_script


    output:
        path ("plots/*.pdf")

    script:
    """
    Rscript $filter_stats_script ${stats_folder}

    """
}