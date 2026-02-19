#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// WE REQUIRE THIS TO BE DONE BEFORE
// download mixer.zip at (https://www.dropbox.com/scl/fo/y5yl2bd5mgplsjwwzsx77/AIFIhSJkzJTFIYhR95TwRVc?rlkey=eydtbzwva5294snzgf6lz0g5f&e=1&st=5g7d4jiv&dl=0)
// mkdir -p /scratch/c.username/mixer_dropbox
// mv ~/mixer.zip /scratch/c.username/mixer_dropbox/

// replace this to you corresponding reference data (e.g. 1,000 genomes - bim/bed/fam)
// The current script is ONLY optimised for 1Kg reference data

// params for reformat sumstats
if( !params.containsKey('beta_col') ) {
    params.beta_col = "BETA"
}
if( !params.containsKey('se_col') ) {
    params.se_col = "SE"
}
if( !params.containsKey('pos_col') ) {
    params.pos_col = "POS"
}

// params for mixer itself 
if( !params.containsKey('mixer') ) {
    params.mixer = "MIXER_PY"
}
if( !params.containsKey('bashbin') ) {
    params.bashbin = "bash"
}
if( !params.containsKey('pybin') ) {
    params.pybin = "python3"
}
if( !params.containsKey('user') ) {
    params.user = ""
}
if( !params.containsKey('do_univ') ) {
    params.do_univ = false
}
if( !params.containsKey('do_biv') ) {
    params.do_biv = false
}
if( !params.containsKey('traits') ) {
    params.traits = ""
}
if( !params.containsKey('phenos') ) {
    params.phenos = ""
}
if( !params.containsKey('sumstats_kind') ) {
    params.sumstats_kind = "postqc"
}
if( !params.containsKey('sumstats_pattern') ) {
    params.sumstats_pattern = ""
}
if( !params.containsKey('pairs') ) {
    params.pairs = ""
}
if( !params.containsKey('rep') ) {
    params.rep = 1
}
if( !params.containsKey('mixer_root') ) {
    params.mixer_root = "${workflow.launchDir}/ref/mixer_dropbox"
}
if( !params.containsKey('ld_ref') ) {
    params.ld_ref = "${params.mixer_root}/reference/ldsc/1000G_EUR_Phase3_plink"
}
if( !params.containsKey('out_base') ) {
    params.out_base = "${workflow.launchDir}/outputs/mixer"
}
if( !params.containsKey('logdir') ) {
    params.logdir = "${workflow.launchDir}/outputs/logs/MiXeR"
}

if( !params.user ) {
    error "Missing required parameter: --user (e.g. --user c.c99999999)"
}

if( params.sumstats_kind != "postqc" && params.sumstats_kind != "mixer_ready" ) {
    error "Unknown --sumstats_kind ${params.sumstats_kind} (use postqc or mixer_ready)"
}

if( !params.traits && !params.phenos ) {
    error "Missing required parameter: --traits traits.csv OR --phenos AD,SCZ,LON"
}

if( !params.sumstats_pattern ) {
    if( params.sumstats_kind == "postqc" ) {
        params.sumstats_pattern = "${workflow.launchDir}/data/{PHENO}/post-qc/{PHENO}.ldsc_ready_neff.tsv"
    }
    if( params.sumstats_kind == "mixer_ready" ) {
        params.sumstats_pattern = "${workflow.launchDir}/data/{PHENO}/mixer_ready/{PHENO}_mixer_ready.tsv.gz"
    }
}

process REFORMAT_SUMSTATS {
    tag "${pheno_id}_sumstats"

    input:
    tuple val(pheno_id),
          path(sumstats),
          val(se_col),
          val(beta_col),
          val(pos_col),
          path(out_dir)

    output:
    tuple val(pheno_id), path("${pheno_id}_mixer_ready.tsv.gz")

    script:
        """
    set -euo pipefail

    if command -v module >/dev/null 2>&1; then
        module purge
        module load singularity
    fi

    if [[ -n "${params.logdir}" ]]; then
        mkdir -p "${params.logdir}"
    fi

    mkdir -p "${workflow.launchDir}/data/${pheno_id}/mixer_ready"
    mkdir -p "${out_dir}"

    ${params.pybin} ${workflow.launchDir}/src/mixer/reformat_sumstats.py \
        --sumstats ${sumstats} \
        --beta_col ${beta_col} \
        --se_col ${se_col} \
        --pos_col ${pos_col} \
        --pheno_id ${pheno_id} \
        --out_dir ${out_dir}

    cp ${out_dir}/${pheno_id}_mixer_ready.tsv.gz ./${pheno_id}_mixer_ready.tsv.gz

    if [[ -n "${params.logdir}" ]]; then
        cp .command.* "${params.logdir}/" || true
    fi
    """
}

process UNIVARIATE_MIXER {
    tag "${pheno_id}_univ_mixer"

    when:
    params.do_univ

    input:
    tuple val(pheno_id),
          path(pheno_sumstats)

    output:
    tuple val(pheno_id),
          path("${pheno_id}.rep${params.rep}.json"),
          path("${pheno_id}.rep${params.rep}.fit.json"),
          path("${pheno_id}.rep${params.rep}.log")

    script:
    """
    set -euo pipefail

    if [[ -n "${params.logdir}" ]]; then
        mkdir -p "${params.logdir}"
    fi

    OUTDIR="${params.out_base}/univariate/${pheno_id}"
    mkdir -p "\${OUTDIR}"

    ${params.bashbin} ${workflow.launchDir}/src/mixer/univ_mixer.sh \
        ${params.user} \
        ${pheno_sumstats} \
        ${pheno_id} \
        --rep ${params.rep} \
        --mixer-root ${params.mixer_root} \
        --mixer-ref ${params.ld_ref} \
        --outdir "\${OUTDIR}"

    cp "\${OUTDIR}/${pheno_id}.rep${params.rep}.json" ./${pheno_id}.rep${params.rep}.json
    cp "\${OUTDIR}/${pheno_id}.rep${params.rep}.fit.json" ./${pheno_id}.rep${params.rep}.fit.json
    cp "\${OUTDIR}/${pheno_id}.rep${params.rep}.log" ./${pheno_id}.rep${params.rep}.log

    if [[ -n "${params.logdir}" ]]; then
        cp .command.* "${params.logdir}/" || true
    fi
    """
}

process BIVARIATE_MIXER {
    tag "${pheno1}_${pheno2}_biv_mixer"

    when:
    params.do_biv

    input:
    tuple val(pheno1),
          val(pheno2),
          path(sum1),
          path(sum2),
          path(u1),
          path(u2)

    output:
    tuple val(pheno1), val(pheno2),
          path("${pheno1}_${pheno2}.fit.rep${params.rep}.json"),
          path("${pheno1}_${pheno2}.apply.rep${params.rep}.json"),
          path("${pheno1}_${pheno2}.fit.rep${params.rep}.log"),
          path("${pheno1}_${pheno2}.apply.rep${params.rep}.log")

    script:
    """
    set -euo pipefail

    if [[ -n "${params.logdir}" ]]; then
        mkdir -p "${params.logdir}"
    fi

    OUTDIR="${params.out_base}/bivariate/${pheno1}_${pheno2}"
    mkdir -p "\${OUTDIR}"

    ${params.bashbin} ${workflow.launchDir}/src/mixer/biv_mixer.sh \
        ${params.user} \
        ${sum1} \
        ${sum2} \
        ${pheno1} \
        ${pheno2} \
        --rep ${params.rep} \
        --mixer-root ${params.mixer_root} \
        --mixer-ref ${params.ld_ref} \
        --univdir ${params.out_base}/univariate \
        --outdir "\${OUTDIR}"

    cp "\${OUTDIR}/${pheno1}_${pheno2}.fit.rep${params.rep}.json" ./${pheno1}_${pheno2}.fit.rep${params.rep}.json
    cp "\${OUTDIR}/${pheno1}_${pheno2}.apply.rep${params.rep}.json" ./${pheno1}_${pheno2}.apply.rep${params.rep}.json
    cp "\${OUTDIR}/${pheno1}_${pheno2}.fit.rep${params.rep}.log" ./${pheno1}_${pheno2}.fit.rep${params.rep}.log
    cp "\${OUTDIR}/${pheno1}_${pheno2}.apply.rep${params.rep}.log" ./${pheno1}_${pheno2}.apply.rep${params.rep}.log

    if [[ -n "${params.logdir}" ]]; then
        cp .command.* "${params.logdir}/" || true
    fi
    """
}

workflow {

        ch_inputs = Channel.empty()

    if( params.traits ) {
        ch_inputs = Channel
            .fromPath(params.traits)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    row.pheno as String,
                    file(row.sumstats as String),
                    params.se_col,
                    params.beta_col,
                    params.pos_col,
                    file("${workflow.launchDir}/data/${row.pheno as String}/mixer_ready")
                )
            }
    }

    if( params.phenos ) {
        phenolist = params.phenos
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }

        ch_inputs = Channel
            .fromList(phenolist)
            .map { pheno ->
                sumstats_path = params.sumstats_pattern.replace("{PHENO}", pheno)
                tuple(
                    pheno as String,
                    file(sumstats_path),
                    params.se_col,
                    params.beta_col,
                    params.pos_col,
                    file("${workflow.launchDir}/data/${pheno}/mixer_ready")
                )
            }
    }

    ch_ready = Channel.empty()

    if( params.sumstats_kind == "postqc" ) {
        ch_ready = REFORMAT_SUMSTATS(ch_inputs)
    }

    if( params.sumstats_kind == "mixer_ready" ) {
        ch_ready = ch_inputs.map { pheno_id, sumstats, se_col, beta_col, pos_col, out_dir ->
            tuple(pheno_id, sumstats)
        }
    }

    ch_univ_jsons = Channel.empty()

    if( params.do_univ ) {
        ch_univ_jsons = UNIVARIATE_MIXER(ch_ready)
    }

    if( params.do_biv ) {

        ch_biv_inputs = Channel.empty()

        if( params.do_univ ) {

            ch_map_univ = ch_univ_jsons.map { pheno_id, json_main, json_fit, logf ->
                tuple(pheno_id, json_main)
            }

            ch_biv_inputs = ch_ready
                .cross(ch_ready)
                .filter { a, b -> a[0].toString() < b[0].toString() }
                .map { a, b -> tuple(a[0], b[0], a[1], b[1]) }
                .combine(ch_map_univ, by: 0)
                .map { ph1, ph2, s1, s2, u1 ->
                    tuple(ph1, ph2, s1, s2, u1)
                }
                .combine(ch_map_univ, by: 1)
                .map { ph1, ph2, s1, s2, u1, u2 ->
                    tuple(ph1, ph2, s1, s2, u1, u2)
                }

        } else {

            error "If --do_biv true, set --do_univ true (bivariate requires univariate jsons)"
        }

        BIVARIATE_MIXER(ch_biv_inputs)
    }
}
