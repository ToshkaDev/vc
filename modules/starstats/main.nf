#!/usr/bin/env nextflow

process STAR_STATS {
    container 'bioconductor/bioconductor_docker:RELEASE_3_21'

    publishDir "results/starstats", mode: 'copy'

    input:
        path logs_files
        path star_stats_script


    output:
        path("*.tsv"), emit: star_stats
        path "plots", optional: true

    script:
    """
    Rscript $star_stats_script ${logs_files.join(" ")}

    """
}