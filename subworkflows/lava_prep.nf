#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { LAVA_GWAS_PREP } from '../modules/lava_prep'
include { MAKE_INFO_FILE } from '../modules/lava_prep'

workflow LAVA_PREP {

    take:
    ch_sumstats
    ch_pairs

    main:
    lava_data_prep = file("${workflow.launchDir}/bin/prep_data.py")

    ch_lava = LAVA_GWAS_PREP(ch_sumstats, lava_data_prep).lava_tsv
    ch_lava_key = ch_lava.map { meta, tsv -> tuple(meta.id, tsv) }

    ch_pairs_with_t1 = ch_pairs
        .join(ch_lava_key, by: 0)
        .map { trait1, trait2, meta, t1_tsv ->
            tuple(trait2, meta, t1_tsv)
        }

    ch_info_in = ch_pairs_with_t1
        .join(ch_lava_key, by: 0)
        .map { trait2, meta, t1_tsv, t2_tsv ->
            tuple(meta, t1_tsv, t2_tsv)
        }

    MAKE_INFO_FILE(ch_info_in)

    emit:
    lava_tsv = ch_lava
    info_tsv = MAKE_INFO_FILE.out.info_tsv
}