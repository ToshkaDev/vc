#!/usr/bin/env nextflow

process FILTER_STATS {
    container 'bioconductor/bioconductor_docker:RELEASE_3_21'

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