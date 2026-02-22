#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { HDL_L } from "../modules/hdl_l.nf"

workflow HDL_L_PAIRS {

    take:
    ch_hdl_pairs

    main:
    ch_hdl = HDL_L(ch_hdl_pairs)

    emit:
    hdl = ch_hdl
}