#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process HDL_L {
  tag "${meta.trait1}_${meta.trait2}_hdl"

  publishDir "${params.outdir}/hdl/${meta.trait1}_${meta.trait2}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), path(trait1_sumstats), path(trait2_sumstats)
  file(hdl_r)
  path(ld_path)
  file(bim_path)

  output:
  tuple val(meta), path("${meta.trait1}_${meta.trait2}_hdl"), emit: hdll_out

  script:
  """
  set -euo pipefail

  RBIN="${params.rbin}"
  if [ -z "\$RBIN" ] || [ "\$RBIN" = "null" ]; then
    RBIN="Rscript"
  fi

  PAIR="${meta.trait1}_${meta.trait2}"
  WORKOUT="\${PAIR}_hdl"
  mkdir -p "\$WORKOUT"

  "\$RBIN" "${hdl_r}" \
    "${trait1_sumstats}" \
    "${trait2_sumstats}" \
    "${meta.trait1}" \
    "${meta.trait2}" \
    "${meta.beta1}" \
    "${meta.se1}" \
    "${meta.beta2}" \
    "${meta.se2}" \
    "${ld_path}" \
    "${bim_path}" \
    "\$WORKOUT" \
    "${params.cov_min}"
  """
}