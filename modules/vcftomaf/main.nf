#!/usr/bin/env nextflow

params.outdir = "results/vcftomaf"

process VCF_TO_MAF {
    container 'ensemblorg/ensembl-vep'

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(hc_pass_vcf), path(hc_pass2_vcf), path(hc_pass2_lcr_vcf), path(hc_pass2_lcr_1kG_vcf), path(hc_pass2_lcr_dbsnp_vcf) 
        path genome
        path vcf2maf_script

    output:
        path "*.tab"
        path "*Log.out"
        path "*Log.final.out", emit: logs
        tuple val(sample_id), path("*.bam"), emit: bam

    script:
    """
    # + convert vcf to maf for genvisr
    # vep: ~/anaconda3/envs/py3.7/bin/ => ensembl-vep-105 => gnomad r2.1.1 >125k exomes
    # vcf2maf.pl => default max_subpop_af > 0.0004 any gnomad population
    # here: --max-subpop-af 0.01 => af > 0.01

    # no filter (rna-edit, lcr sites)
    tmp=${hc_pass_vcf/.gz/}
    gunzip -c ${hc_pass_vcf} > $tmp
    perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $tmp --output-maf $tmp.maf --ref-fasta $genome --tumor-id ${sample_id} --max-subpop-af 0.01 --ncbi-build GRCh38 --vep-path /opt/vep/src/ensembl-vep
    pigz -p 4 $tmp.maf
    pigz -p 4 ${tmp/vcf/vep.vcf} 
    rm $tmp

    # with lcr, without rna-edit sites
    tmp_pass2=${hc_pass2_vcf/.gz/}
    gunzip -c ${hc_pass2_vcf} > ${tmp_pass2}
    perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf ${tmp_pass2} --output-maf ${tmp_pass2}.maf --ref-fasta $genome --tumor-id ${sample_id} --ncbi-build GRCh38 --max-subpop-af 0.01 --vep-path /opt/vep/src/ensembl-vep
    pigz -p 4 ${tmp_pass2}.maf
    pigz -p 4 ${tmp_pass2/vcf/vep.vcf} 
    rm ${tmp_pass2}

    # without lcr, without rna-edit sites
    tmp_pass2_lcr=${hc_pass2_lcr_vcf/.gz/}
    gunzip -c ${hc_pass2_lcr_vcf} > ${tmp_pass2_lcr}
    perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf ${tmp_pass2_lcr} --output-maf ${tmp_pass2_lcr}.maf --ref-fasta $genome --tumor-id ${sample_id} --ncbi-build GRCh38 --max-subpop-af 0.01 --vep-path /opt/vep/src/ensembl-vep
    pigz -p 4 ${tmp_pass2_lcr}.maf
    pigz -p 4 ${tmp_pass2_lcr/vcf/vep.vcf} 
    rm ${tmp_pass2_lcr}

    # without lcr, without rna-edit sites, gatk 1000 genomes
    tmp_pass2_lcr_1kG=${hc_pass2_lcr_1kG_vcf/.gz/}
    gunzip -c ${hc_pass2_lcr_1kG_vcf} > ${tmp_pass2_lcr_1kG}
    perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf ${tmp_pass2_lcr_1kG} --output-maf ${tmp_pass2_lcr_1kG}.maf --ref-fasta $genome --tumor-id ${sample_id} --ncbi-build GRCh38 --max-subpop-af 0.01 --vep-path /opt/vep/src/ensembl-vep
    pigz -p 4 ${tmp_pass2_lcr_1kG}.maf
    pigz -p 4 ${tmp_pass2_lcr_1kG/vcf/vep.vcf} 
    rm ${tmp_pass2_lcr_1kG}

    # without lcr, without rna-edit sites, nih dbsnp
    tmp_pass2_lcr_dbsnp=${hc_pass2_lcr_dbsnp_vcf/.gz/}
    gunzip -c ${hc_pass2_lcr_dbsnp_vcf} > ${tmp_pass2_lcr_dbsnp}
    perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf ${tmp_pass2_lcr_dbsnp} --output-maf ${tmp_pass2_lcr_dbsnp}.maf --ref-fasta $genome --tumor-id ${sample_id} --ncbi-build GRCh38 --max-subpop-af 0.01 --vep-path /opt/vep/src/ensembl-vep
    pigz -p 4 ${tmp_pass2_lcr_dbsnp}.maf
    pigz -p 4 ${tmp_pass2_lcr_dbsnp/vcf/vep.vcf} 
    rm ${tmp_pass2_lcr_dbsnp}

    """
}

process COMPRESS_MAF {
    container 'biocontainers/pigz:v2.3.4_cv1'

    input:
        tuple val(sample_id), path(maf), path(vep), path(maf_pass2), path(vep_pass2), path(maf_pass2_lcr), path(vep_pass2_lcr), 
        path(maf_pass2_lcr_1kG), path(vep_pass2_lcr_1kG), path(maf_pass2_lcr_dbsnp), path(vep_pass2_lcr_dbsnp)

    output:
        tuple val(sample_id), path("*.gz")
        // tuple val(sample_id), path("${maf.name}.gz"), path("${vep.name}.gz"), path("${maf_pass2.name}.gz"), path("${vep_pass2.name}.gz"),
        // path("${maf_pass2_lcr.name}.gz"), path("${vep_pass2_lcr.name}.gz"), path("${maf_pass2_lcr_1kG.name}.gz"), path("${vep_pass2_lcr_1kG.name}.gz"),
        // path("${maf_pass2_lcr_dbsnp.name}.gz"), path("${vep_pass2_lcr_dbsnp.name}.gz")


    script:
    """
    pigz -p 4 $maf
    pigz -p 4 $vep
    pigz -p 4 $maf_pass2
    pigz -p 4 $vep_pass2
    pigz -p 4 $maf_pass2_lcr
    pigz -p 4 $vep_pass2_lcr
    pigz -p 4 $maf_file
    pigz -p 4 $vep_vcf_file
    """
}
