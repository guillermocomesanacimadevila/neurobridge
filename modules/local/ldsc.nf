#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LDSC {
  tag "${meta.trait1}_${meta.trait2}"

  publishDir "${params.outdir}/ldsc/${meta.trait1}_${meta.trait2}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), path(trait1_sumstats), path(trait2_sumstats)
  file(ldsc_r)
  file(hm3_snplist)
  path(ld_chr_dir)
  path(wld_dir)

  output:
  tuple val(meta), path("${meta.trait1}_${meta.trait2}_ldsc"), emit: ldsc_outdir

  script:
  """
  set -euo pipefail

  RBIN="${params.rbin}"
  if [ -z "\$RBIN" ] || [ "\$RBIN" = "null" ]; then
    RBIN="Rscript"
  fi

  PAIR="${meta.trait1}_${meta.trait2}"
  WORKOUT="\${PAIR}_ldsc"
  mkdir -p "\$WORKOUT"

  "\$RBIN" "${ldsc_r}" \
    "${trait1_sumstats}" \
    "${trait2_sumstats}" \
    "${meta.trait1}" \
    "${meta.trait2}" \
    "${meta.cases1}" \
    "${meta.controls1}" \
    "${meta.cases2}" \
    "${meta.controls2}" \
    "${meta.pop_prev1}" \
    "${meta.pop_prev2}" \
    "\$WORKOUT" \
    "${hm3_snplist}" \
    "${ld_chr_dir}" \
    "${wld_dir}"
  """
}