#!/usr/bin/env nextflow

params.outdir = "results/vcftomaf"

process VCF_TO_MAF {
    container 'community.wave.seqera.io/library/ensembl-vep_vcf2maf:6836d39d9a125f9f'

    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), path(hc_pass_vcf), path(_), path(hc_pass2_vcf), path(hc_pass2_lcr_vcf), path(hc_pass2_lcr_1kG_vcf), path(hc_pass2_lcr_dbsnp_vcf) 
        path genome
        path vep_data_cache

    output:
        path "*"

    script:
    """
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

    vcf_to_maf "${hc_pass_vcf}"
    """
}
