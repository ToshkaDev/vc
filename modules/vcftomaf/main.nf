#!/usr/bin/env nextflow

params.outdir = "results/vcftomaf"

process VCF_TO_MAF {
    container 'community.wave.seqera.io/library/ensembl-vep_vcf2maf:6836d39d9a125f9f'

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(hc_pass_vcf), path(_), path(hc_pass2_vcf), path(hc_pass2_lcr_vcf), 
            path(hc_pass2_lcr_1kG_vcf), path(hc_pass2_lcr_dbsnp_vcf) 
        path genome
        path vep_data_cache

    output:
        tuple val(sample_id), path("*.hc.pass.vcf.maf"), path("*.hc.pass.vep.vcf"), path("*.hc.pass2.vcf.maf"), 
            path("*.hc.pass2.vep.vcf"), path("*.hc.pass2.lcr.vcf.maf"), path("*.hc.pass2.lcr.vep.vcf"), 
            path("*.hc.pass2.lcr.1kG.vcf.maf"), path("*.hc.pass2.lcr.1kG.vep.vcf"), path("*.hc.pass2.lcr.dbsnp.vcf.maf"), 
            path("*.hc.pass2.lcr.dbsnp.vep.vcf"), emit: all_maf_vep

    script:
    """
    # + convert vcf to maf for genvisr
    # vep: ~/anaconda3/envs/py3.7/bin/ => ensembl-vep-105 => gnomad r2.1.1 >125k exomes
    # vcf2maf.pl => default max_subpop_af > 0.0004 any gnomad population
    # here: --max-subpop-af 0.01 => af > 0.01
    
    vcf_to_maf() {
        vcf=\$1
        tmpdir=/tmp/tmp_dir
        mkdir -p \$tmpdir

        #tmp=\${vcf%.*}
        tmp_basename=\$(basename "\${vcf%.gz}")
        tmp="\$tmpdir/\$tmp_basename"

        gunzip -c \${vcf} > \$tmp
        perl /opt/conda/bin/vcf2maf.pl \
        --input-vcf \$tmp \
        --output-maf \${tmp}.maf \
        --ref-fasta $genome \
        --tumor-id ${sample_id} \
        --max-subpop-af 0.01 \
        --ncbi-build GRCh38 \
        --vep-path /opt/conda/bin \
        --vep-data ${vep_data_cache} \
        --tmp-dir \$tmpdir

        # Move the intermediate VEP VCF and MAF out before cleanup
        mv \$tmpdir/*.vep.vcf .
        mv \$tmpdir/*.maf .
        rm \$tmp       
    }

    #no filter (rna-edit, lcr sites)
    vcf_to_maf "${hc_pass_vcf}"
    
    #with lcr, without rna-edit sites
    vcf_to_maf "${hc_pass2_vcf}"
    
    #without lcr, without rna-edit sites
    vcf_to_maf "${hc_pass2_lcr_vcf}"
    
    #without lcr, without rna-edit sites, gatk 1000 genomes
    vcf_to_maf "${hc_pass2_lcr_1kG_vcf}"

    #without lcr, without rna-edit sites, nih dbsnp
    vcf_to_maf "${hc_pass2_lcr_dbsnp_vcf}"
    """
}

process COMPRESS_MAF_VEP {
    container 'quay.io/biocontainers/pigz:2.8'

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(pass_vcf_maf), path(pass_vep), path(pass2_vcf_maf), path(pass2_vep), 
            path(pass2_lcr_vcf_maf), path(pass2_lcr_vep), path(pass2_lcr_1kG_vcf_maf), 
        path(pass2_lcr_1kG_vep), path(pass2_lcr_dbsnp_vcf_maf), path(pass2_lcr_dbsnp_vep)

    output:
        // tuple val(sample_id), path("${pass_vcf_maf.name}.gz"), path("${pass_vep.name}.gz"), emit: pass_vcf_gz
        // tuple val(sample_id), path("${pass2_vcf_maf.name}.gz"), path("${pass2_vep.name}.gz"), emit: pass2_vcf_gz
        // tuple val(sample_id), path("${pass2_lcr_vcf_maf.name}.gz"), path("${pass2_lcr_vep.name}.gz"), emit: pass2_lcr_vcf_gz
        // tuple val(sample_id), path("${pass2_lcr_1kG_vcf_maf.name}.gz"), path("${pass2_lcr_1kG_vep.name}.gz"), emit: pass2_lcr_1kG_vcf_gz
        // tuple val(sample_id), path("${pass2_lcr_dbsnp_vcf_maf.name}.gz"), path("${pass2_lcr_dbsnp_vep.name}.gz"), emit: pass2_lcr_dbsnp_vcf_gz
        tuple val(sample_id), path("${pass_vcf_maf.name}.gz"), path("${pass2_vcf_maf.name}.gz"), path("${pass2_lcr_vcf_maf.name}.gz"), 
            path("${pass2_lcr_1kG_vcf_maf.name}.gz"), path("${pass2_lcr_dbsnp_vcf_maf.name}.gz"), emit: compressed_maff
        path("*.vep.vcf.gz")

    script:
    """
    compress() {
        for f in "\$@"; do
            realfile=\$(readlink -f "\$f")
            cp "\$realfile" "./tmp_\${f}"
            pigz -p 4 "./tmp_\${f}"
            mv "./tmp_\${f}.gz" "\${f}.gz"
            rm "\$f"
        done
    }
    compress $pass_vcf_maf $pass_vep
    compress $pass2_vcf_maf $pass2_vep
    compress $pass2_lcr_vcf_maf $pass2_lcr_vep
    compress $pass2_lcr_1kG_vcf_maf $pass2_lcr_1kG_vep
    compress $pass2_lcr_dbsnp_vcf_maf $pass2_lcr_dbsnp_vep
    """
}
