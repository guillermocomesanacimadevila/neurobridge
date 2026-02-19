#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// WE ASSUME PRE-PROCESS LDSC data
// SO RUN ldsc.nf first and this script USES that data as input
// reference files - zenodo 

params.cov_min = 0.90 // coverage
params.rbin = "Rscript"
params.ld_path = "/Users/c24102394/ref/HDL-L_ref/LD.path"
params.bim_path = "/Users/c24102394/ref/HDL-L_ref/bimfile"

process HDL_L {
    tag "${pheno1_prefix}_${pheno2_prefix}_hdl_l"

    input:
    tuple  path(gwas1),
           path(gwas2),
           val(pheno1_prefix),
           val(pheno2_prefix),
           val(beta1_col),
           val(se1_col),
           val(beta2_col),
           val(se2_col),
           val(out_dir)

    output:
    path "HDLL_out", emit: hdll_out

    script:
    """
    set -euo pipefail

    mkdir -p HDLL_out

    ${params.rbin} ${workflow.launchDir}/src/HDL-L/hdl_l.R \
        ${gwas1} \
        ${gwas2} \
        ${pheno1_prefix} \
        ${pheno2_prefix} \
        ${beta1_col} \
        ${se1_col} \
        ${beta2_col} \
        ${se2_col} \
        ${params.ld_path} \
        ${params.bim_path} \
        HDLL_out \
        ${params.cov_min}

    mkdir -p ${workflow.launchDir}/outputs/HDL-L/${pheno1_prefix}-${pheno2_prefix}
    cp -r HDLL_out/* ${workflow.launchDir}/outputs/HDL-L/${pheno1_prefix}-${pheno2_prefix}/

    mkdir -p ${workflow.launchDir}/logs/HDL_L
    cp .command.* ${workflow.launchDir}/logs/HDL_L/ || true
    """
}

workflow {

    traits = ["AD", "SCZ", "LON"]

    trait_pairs = []
    for (i in 0..<(traits.size()-1)) {
        for (j in (i+1)..<(traits.size())) {
            trait_pairs << [traits[i], traits[j]]
        }
    }

    pairs = Channel
        .fromList(trait_pairs)
        .map { pair ->
            def a = pair[0]
            def b = pair[1]
            def gwas1_path = "${workflow.launchDir}/data/Main/${a}/post-qc/${a}.ldsc_ready_neff.tsv"
            def gwas2_path = "${workflow.launchDir}/data/Main/${b}/post-qc/${b}.ldsc_ready_neff.tsv"
            def out_dir = "${workflow.launchDir}/outputs/HDL-L/${a}-${b}"

            tuple(
                file(gwas1_path),
                file(gwas2_path),
                a,
                b,
                "BETA",
                "SE",
                "BETA",
                "SE",
                out_dir
            )
        }

    HDL_L(pairs)
}
