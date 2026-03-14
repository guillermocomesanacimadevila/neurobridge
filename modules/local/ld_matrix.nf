#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GEN_LD_MATRIX {

    tag "${meta.trait1}_${meta.trait2}_gen_ld_matrix"

    input:
    tuple val(meta), path(loci_dir)
    path gen_ld_bash
    path ref_chr_files, stageAs: 'refpanel/*'
    val  ref_chr_base

    output:
    tuple val(meta), path("defined_loci"), emit: ld_ready

    script:
    """
    set -euo pipefail

    BASHBIN="bash"

    "\$BASHBIN" "${gen_ld_bash}" \\
        "${meta.trait1}" \\
        "${meta.trait2}" \\
        "${loci_dir}" \\
        "refpanel/${ref_chr_base}"
    """
}