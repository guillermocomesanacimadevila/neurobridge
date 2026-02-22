#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { LDSC } from "../modules/ldsc"

workflow LDSC_PAIRWISE {
  take:
  ch_pairs

  main:
  ch_ldsc = LDSC(ch_pairs)

  emit:
  ldsc = ch_ldsc
}