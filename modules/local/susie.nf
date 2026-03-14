#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SUSIE_OVERLAP_MAP {

    tag "${meta.pair}_${meta.locus}_L${L}"
    publishDir "${params.outdir}/SuSiE/${meta.pair}/${meta.locus}/L${L}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(gwas1), path(gwas2), path(ld_matrix)
    path susie_r
    path overlap_py
    path map_gene_py
    path gtf
    val L

    output:
    tuple val(meta), path("${meta.trait1}_${meta.trait2}_${meta.coords}.tsv"), emit: overlap
    tuple val(meta), path("mapped_${meta.trait1}_${meta.trait2}_${meta.coords}.tsv"), emit: mapped
    tuple val(meta), path("cs95_${meta.trait1}_L${L}.tsv"), emit: cs95_trait1
    tuple val(meta), path("cs95_${meta.trait2}_L${L}.tsv"), emit: cs95_trait2
    tuple val(meta), path("pip_${meta.trait1}_L${L}.tsv"), emit: pip_trait1
    tuple val(meta), path("pip_${meta.trait2}_L${L}.tsv"), emit: pip_trait2
    tuple val(meta), path("res_${meta.trait1}_L${L}.rds"), emit: rds_trait1
    tuple val(meta), path("res_${meta.trait2}_L${L}.rds"), emit: rds_trait2

    script:
    """
    set -euo pipefail

    RBIN="${params.rbin}"
    if [ -z "\$RBIN" ] || [ "\$RBIN" = "null" ]; then
        RBIN="Rscript"
    fi

    PYBIN="${params.pybin}"
    if [ -z "\$PYBIN" ] || [ "\$PYBIN" = "null" ]; then
        PYBIN="python3"
    fi

    # trait1 SuSiE
    "\$RBIN" "${susie_r}" \\
        "${meta.trait1}" \\
        "${gwas1}" \\
        "${L}" \\
        "${ld_matrix}" \\
        "${meta.n1}" \\
        "${params.max_iter}"

    # trait2 SuSiE
    "\$RBIN" "${susie_r}" \\
        "${meta.trait2}" \\
        "${gwas2}" \\
        "${L}" \\
        "${ld_matrix}" \\
        "${meta.n2}" \\
        "${params.max_iter}"

    # overlap cs95 sets
    "\$PYBIN" "${overlap_py}" \\
        --pheno1_id "${meta.trait1}" \\
        --pheno2_id "${meta.trait2}" \\
        --pheno1_path "cs_${meta.trait1}_L${L}.tsv" \\
        --pheno2_path "cs_${meta.trait2}_L${L}.tsv" \\
        --out_path "${meta.trait1}_${meta.trait2}_${meta.coords}.tsv"

    # map overlap SNPs to closest gene
    "\$PYBIN" "${map_gene_py}" \\
        --input "${meta.trait1}_${meta.trait2}_${meta.coords}.tsv" \\
        --gtf "${gtf}" \\
        --out "mapped_${meta.trait1}_${meta.trait2}_${meta.coords}.tsv"
    """
}