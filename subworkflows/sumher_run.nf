#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TAGGING } from '../modules/tagging'
include { H2 }      from '../modules/h2'
include { RG }      from '../modules/rg'

workflow SUMHER_RUN {

    take:
    ch_sumstats
    ch_pairs
    plink_dir
    calc_p
    ldak_bin

    main:
    def tagfile

    if( params.do_tagging ) {
        tagfile = TAGGING(ldak_bin, plink_dir).tagfiles.map { tagging, taglist -> tagging }.first()
    } else {
        tagfile = Channel.fromPath("${params.outdir}/sumher/tagging/${params.tag_prefix}.tagging", checkIfExists: true).first()
    }

    ch_h2 = H2(ch_sumstats, tagfile, ldak_bin)
    ch_summaries = ch_h2.sums.map { meta, sums -> tuple(meta.id, sums) }

    ch_with_t1 = ch_pairs
        .join(ch_summaries, by: 0)
        .map { trait1, trait2, meta, s1 -> tuple(trait2, meta, s1) }

    ch_rg_in = ch_with_t1
        .join(ch_summaries, by: 0)
        .map { trait2, meta, f1, f2 -> tuple(meta, f1, f2) }

    RG(ch_rg_in, tagfile, ldak_bin, calc_p)

    emit:
    h2_summaries = ch_h2.sums
}