#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/neurobridge -> LAVA REMINDERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Bloke needs to download reference data from -> dropbox
    * Bloke needs to run LDSC workflow/ before running LAVA 
    * Run ./workflows/neurobridge/main_ldsc.nf

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process LAVA {
    tag "${meta.trait1}_${meta.trait2}_lava_run"
    publishDir "${params.outdir}/lava/res/${meta.trait1}_${meta.trait2}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(t1_tsv), path(t2_tsv), path(info_tsv), path(overlap_csv)
    path(lava_r)
    path(lava_ref_dir)
    path(loci_file)

    output:
    tuple val(meta), path("LAVA_local_rg_bivariate.tsv"), emit: lava_out
    path("LAVA_local_rg_partial.tsv"), emit: lava_partial

    script:
    """
    set -euo pipefail
    REF_PREFIX="${lava_ref_dir}/lava-ukb-v1.1"
    Rscript "${lava_r}" \\
      "\$REF_PREFIX" \\
      "${loci_file}" \\
      "${info_tsv}" \\
      "${overlap_csv}" \\
      "." \\
      "${meta.trait1}" "${meta.cases1}" "${meta.controls1}" "${meta.pop_prev1}" "${t1_tsv}" \\
      "${meta.trait2}" "${meta.cases2}" "${meta.controls2}" "${meta.pop_prev2}" "${t2_tsv}"
    """
}