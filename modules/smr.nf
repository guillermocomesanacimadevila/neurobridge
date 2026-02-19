#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.pybin = "python3"
params.ad_sumstats = "${workflow.launchDir}/data/AD/post-qc/AD.ldsc_ready_neff.tsv"
params.scz_sumstats = "${workflow.launchDir}/data/SCZ/post-qc/SCZ.ldsc_ready_neff.tsv"
params.lon_sumstats = "${workflow.launchDir}/data/LON/post-qc/LON.ldsc_ready_neff.tsv"
params.snp_col = "SNP"
params.a1_col = "A1"
params.a2_col = "A2"
params.freq_col = "FRQ"
params.beta_col = "BETA"
params.se_col = "SE"
params.p_col = "P"
params.n_col = "N"
params.ld_ref = "${workflow.launchDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL"
params.eqtl_ref = "${workflow.launchDir}/ref/eqtl/GTEx_v8" // GTEx v8 tissues only 

process REFORMAT_SUMSTATS {
    tag "${pheno_id}_sMR_prep"

    publishDir "${workflow.launchDir}/outputs/sMR/smr_inputs", mode: 'copy', overwrite: true

    input:
    tuple val(pheno_id),
          path(pheno_sumstats)

    output:
    tuple val(pheno_id),
          path("${pheno_id}.smr.ma")

    script:
    """
    set -euo pipefail
    mkdir -p ${workflow.launchDir}/logs/sMR
    ${params.pybin} ${workflow.launchDir}/src/sMR/format_gwas.py ${pheno_id} \\
        --infile ${pheno_sumstats} \\
        --out_dir . \\
        --snp_col ${params.snp_col} \\
        --a1_col ${params.a1_col} \\
        --a2_col ${params.a2_col} \\
        --freq_col ${params.freq_col} \\
        --beta_col ${params.beta_col} \\
        --se_col ${params.se_col} \\
        --p_col ${params.p_col} \\
        --n_col ${params.n_col}

    cp .command.* ${workflow.launchDir}/logs/sMR/ || true
    """
}

process RUN_ALL_SMR {

    tag "RUN_ALL_SMR"

    input:
    tuple path(ad_ma),
          path(scz_ma),
          path(lon_ma)

    output:
    path("RUN_ALL_SMR.done")

    script:
    """
    set -euo pipefail 
    source ~/.zshrc
    which ${params.smr_bin}
    ${params.smr_bin} --version
    mkdir -p ${workflow.launchDir}/logs/sMR
    chmod +x ${workflow.launchDir}/src/sMR/smr.sh
    chmod +x ${workflow.launchDir}/src/sMR/run_all_smr.sh
    bash ${workflow.launchDir}/src/sMR/run_all_smr.sh \\
        ${ad_ma} \\
        ${scz_ma} \\
        ${lon_ma} \\
        ${params.ld_ref} \\
        ${params.eqtl_ref}
    touch RUN_ALL_SMR.done
    cp .command.* ${workflow.launchDir}/logs/sMR/ || true
    """
}

workflow {
    Channel.of(
        tuple("AD", file(params.ad_sumstats)),
        tuple("SCZ", file(params.scz_sumstats)),
        tuple("LON", file(params.lon_sumstats))
    ) | REFORMAT_SUMSTATS | set { gwas_ch }

    ad_ma  = gwas_ch.filter { it[0] == "AD"  }.map { it[1] }
    scz_ma = gwas_ch.filter { it[0] == "SCZ" }.map { it[1] }
    lon_ma = gwas_ch.filter { it[0] == "LON" }.map { it[1] }

    ad_ma
        .combine(scz_ma)
        .combine(lon_ma)
        .map { ad, scz, lon -> tuple(ad, scz, lon) }
        | RUN_ALL_SMR
}
