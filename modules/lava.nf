#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.pybin = "python3"
params.rbin = "Rscript"

process data_prep {
    tag "${trait}_prep4_lava"
    input:
    tuple val(trait), val(input_pheno)
    output:
    path("${trait}.tsv"), emit: lava_tsv
    script:
    """
    set -e
    mkdir -p ${workflow.launchDir}/data/Main/${trait}/lava-ready
    ${params.pybin} ${workflow.launchDir}/src/lava/prep_data.py \\
        --input ${workflow.launchDir}/data/Main/${trait}/post-qc/${input_pheno} \\
        --output ${workflow.launchDir}/data/Main/${trait}/lava-ready/${trait}.tsv
    cp ${workflow.launchDir}/data/Main/${trait}/lava-ready/${trait}.tsv ${trait}.tsv
    """
}

process MAKE_INFO3 {
    tag "AD_LON_SCZ_info"
    input:
    tuple val(ph1_name), val(ph1_cases), val(ph1_ctrls), val(ph1_prev), val(ph1_relpath),
          val(ph2_name), val(ph2_cases), val(ph2_ctrls), val(ph2_prev), val(ph2_relpath),
          val(ph3_name), val(ph3_cases), val(ph3_ctrls), val(ph3_prev), val(ph3_relpath)
    output:
    path("input.info.txt")
    script:
    """
    set -e
    cat << EOF > input.info.txt
phenotype\tcases\tcontrols\tprevalence\tfilename
${ph1_name}\t${ph1_cases}\t${ph1_ctrls}\t${ph1_prev}\t${ph1_relpath}
${ph2_name}\t${ph2_cases}\t${ph2_ctrls}\t${ph2_prev}\t${ph2_relpath}
${ph3_name}\t${ph3_cases}\t${ph3_ctrls}\t${ph3_prev}\t${ph3_relpath}
EOF
    """
}

process LAVA3 {
    tag "AD_LON_SCZ_LAVA"
    cpus 1
    publishDir "${workflow.launchDir}/outputs/LAVA/AD_LON_SCZ", mode: 'copy', overwrite: true
    input:
    val(lava_ref)
    val(loci_file)
    path(info_tsv)
    val(ov12_file)
    val(ov13_file)
    val(ov23_file)
    path(ad_tsv)
    path(lon_tsv)
    path(scz_tsv)
    output:
    tuple path("LAVA_local_rg_bivariate.tsv"), path("LAVA_local_rg_partial.tsv"), emit: lava3_out
    script:
    """
    set -e
    mkdir -p ${workflow.launchDir}/logs/lava

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export VECLIB_MAXIMUM_THREADS=1
    export DATA_TABLE_NUM_THREADS=1

    ${params.rbin} ${workflow.launchDir}/src/lava/lava.R \\
        ${lava_ref} \\
        ${loci_file} \\
        ${info_tsv} \\
        ${ov12_file} \\
        ${ov13_file} \\
        ${ov23_file} \\
        . \\
        AD 35274 59163 0.07 ${ad_tsv} \\
        LON 354854 354855 0.10 ${lon_tsv} \\
        SCZ 67390 94015 0.01 ${scz_tsv}

    cp .command.* ${workflow.launchDir}/logs/lava/ || true
    """
}

process SPLIT_LAVA_PAIR {
    tag "${pair}"
    cpus 1
    publishDir { pair_out }, mode: 'copy', overwrite: true
    input:
    tuple val(pair), val(pair_out), path(bivar_tsv), path(partial_tsv)
    output:
    path("LAVA_local_rg_bivariate.tsv")
    path("LAVA_local_rg_partial.tsv")
    script:
    """
    set -e
    ${params.pybin} - << 'PY'
import pandas as pd

bivar = pd.read_csv("${bivar_tsv}", sep="\\t")
partial = pd.read_csv("${partial_tsv}", sep="\\t")

if "pair" in bivar.columns:
    bivar = bivar[bivar["pair"] == "${pair}"]
else:
    bivar = bivar.iloc[0:0]

if "pair" in partial.columns:
    partial = partial[partial["pair"] == "${pair}"]
else:
    partial = partial.iloc[0:0]

bivar.to_csv("LAVA_local_rg_bivariate.tsv", sep="\\t", index=False)
partial.to_csv("LAVA_local_rg_partial.tsv", sep="\\t", index=False)
PY
    """
}

workflow {
    trait_info_prep = Channel.of(
        tuple("AD",  "AD.ldsc_ready_neff.tsv"),
        tuple("SCZ", "SCZ.ldsc_ready_neff.tsv"),
        tuple("LON", "LON.ldsc_ready_neff.tsv")
    )

    prep_out = data_prep(trait_info_prep)

    ad_file = prep_out.filter { it.getName() == 'AD.tsv' }.first()
    scz_file = prep_out.filter { it.getName() == 'SCZ.tsv' }.first()
    lon_file = prep_out.filter { it.getName() == 'LON.tsv' }.first()

    lava_ref = "/Users/c24102394/ref/lava/lava_ref/lava-ukb-v1.1"
    loci_file = "/Users/c24102394/ref/lava/hdll_blocks.coords.loci"

    ov_ad_lon = "${workflow.launchDir}/outputs/ldsc/AD_LON/overlap_corr_for_LAVA_AD_LON.csv"
    ov_ad_scz  = "${workflow.launchDir}/outputs/ldsc/AD_SCZ/overlap_corr_for_LAVA_AD_SCZ.csv"
    ov_scz_lon = "${workflow.launchDir}/outputs/ldsc/SCZ_LON/overlap_corr_for_LAVA_SCZ_LON.csv"

    info_in = Channel.of(
        tuple(
            "AD", 35274, 59163, 0.07, "AD.tsv",
            "LON", 354854, 354855, 0.10, "LON.tsv",
            "SCZ", 67390, 94015, 0.01, "SCZ.tsv"
        )
    )

    info_tsv = MAKE_INFO3(info_in)

    lava3 = LAVA3(
        lava_ref,
        loci_file,
        info_tsv,
        ov_ad_lon,
        ov_ad_scz,
        ov_scz_lon,
        ad_file,
        lon_file,
        scz_file
    )

    pair_jobs = Channel.of(
        tuple("AD_SCZ", "${workflow.launchDir}/outputs/LAVA/AD_SCZ"),
        tuple("AD_LON", "${workflow.launchDir}/outputs/LAVA/AD_LON"),
        tuple("SCZ_LON", "${workflow.launchDir}/outputs/LAVA/SCZ_LON")
    )

    split_args = pair_jobs.combine(lava3).map { pj, res ->
        def (pair, pair_out) = pj
        def (bivar_tsv, partial_tsv) = res
        tuple(pair, pair_out, bivar_tsv, partial_tsv)
    }

    SPLIT_LAVA_PAIR(split_args)
}
