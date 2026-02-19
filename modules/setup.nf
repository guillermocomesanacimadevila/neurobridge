#!/usr/bin/env nextlflow
nextflow.enable.dsl=2

// check if py installed
// check if R installed
// check if docker installed
// open -a docker automatically 
// check if conda installed
// check if docker image built or not 

process PREFLIGHT_CHECKS {
  tag "preflight"

  output:
  path "preflight.ok", emit: preflight_ok

  script:
  """
  set -euo pipefail
  command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is not installed or not in PATH" >&2; exit 1; }
  command -v R       >/dev/null 2>&1 || { echo "ERROR: R is not installed or not in PATH" >&2; exit 1; }
  command -v conda   >/dev/null 2>&1 || { echo "ERROR: Conda is not installed or not in PATH" >&2; exit 1; }
  command -v docker  >/dev/null 2>&1 || { echo "ERROR: Docker is not installed or not in PATH" >&2; exit 1; }

  touch preflight.ok
  """
}

workflow SETUP {
  main:
    PREFLIGHT_CHECKS()
}

// include this in main.nf ---->
// include { SETUP } from './modules/setup.nf'

// workflow {
//   SETUP()
// }
