#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { LAVA } from '../modules/lava'

workflow LAVA_RUN {

  take:
  ch_lava_in
  lava_r
  lava_ref_dir
  loci_file

  main:
  LAVA(
    ch_lava_in,
    lava_r,
    lava_ref_dir,
    loci_file
  )

  emit:
  lava_out = LAVA.out.lava_out
  lava_partial = LAVA.out.lava_partial
}