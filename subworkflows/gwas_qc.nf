#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS }  from '../modules/qc_gwas'
include { ADD_NEFF } from '../modules/add_neff'

workflow GWAS_QC {
  take:
  ch_in

  main:
  ch_qc = QC_GWAS(ch_in).ldsc_ready
  ch_neff = ADD_NEFF(ch_qc).ldsc_neff

  emit:
  ldsc_ready = ch_qc
  ldsc_neff = ch_neff
}