#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.out_base="${workflow.launchDir}/outputs/fuma_prep"
params.post_base="${workflow.launchDir}/outputs/fuma_post"
params.clump_base="${workflow.launchDir}/outputs/clumping"

process PREP_4_FUMA {
  tag "${set_id}"
  publishDir "${params.out_base}", mode:'copy'

  input:
  tuple val(set_id), val(clump_dir), val(traits)

  output:
  path("${set_id}")

  script:
  def traits_str=traits.join(' ')
  """
  set -euo pipefail
  mkdir -p ${set_id}/loci ${set_id}/meta
  echo -e "set_id\\ttrait\\tlocus_id\\tlocus_dir\\tfile_gz" > ${set_id}/meta/manifest.tsv
  for t in ${traits_str}; do
    mkdir -p ${set_id}/loci/${'$'}t
    for d in ${clump_dir}/locus_*; do
      [ -d "${'$'}d" ] || continue
      in="${'$'}d/gwas_${'$'}t.tsv"
      [ -f "${'$'}in" ] || continue
      locus_dir="${'$'}(basename "${'$'}d")"
      locus_id="${'$'}(echo "${'$'}locus_dir" | cut -d'_' -f2)"
      out="${set_id}/loci/${'$'}t/${'$'}{locus_dir}_gwas_${'$'}t.tsv.gz"
      python3 - "${'$'}in" "${'$'}out" <<'PY'
import sys,gzip
infile=sys.argv[1]
outfile=sys.argv[2]
NL=chr(10)
with open(infile) as f, gzip.open(outfile,'wt') as out:
  header=f.readline().rstrip().split('\t')
  out.write('\t'.join(header)+NL)
  idx={name:i for i,name in enumerate(header)}
  for line in f:
    row=line.rstrip().split('\t')
    for col in ('N','CHR','POS'):
      if col in idx and row[idx[col]]!='':
        row[idx[col]]=str(int(float(row[idx[col]])))
    out.write('\t'.join(row)+NL)
PY

      echo -e "${set_id}\\t${'$'}t\\t${'$'}locus_id\\t${'$'}locus_dir\\t${'$'}out" >> ${set_id}/meta/manifest.tsv
    done
  done

  [ -f ${clump_dir}/lead_snps.tsv ] && cp ${clump_dir}/lead_snps.tsv ${set_id}/meta/
  [ -f ${clump_dir}/locus_coords.tsv ] && cp ${clump_dir}/locus_coords.tsv ${set_id}/meta/
  """
}

process FUMA_POST_DIRS {
  tag "${trait}"
  publishDir "${params.post_base}", mode:'copy'

  input:
  tuple val(trait), val(clump_dir)

  output:
  path("${trait}")

  script:
  """
  set -euo pipefail
  mkdir -p ${trait}
  for d in ${clump_dir}/locus_*; do
    [ -d "${'$'}d" ] || continue
    locus_dir="${'$'}(basename "${'$'}d")"
    locus_id="${'$'}(echo "${'$'}locus_dir" | cut -d'_' -f2)"
    mkdir -p ${trait}/locus_${'$'}{locus_id}
  done
  """
}

workflow {

  Channel.of(
    tuple("AD-SCZ","${params.clump_base}/AD_SCZ/AD-SCZ",["AD","SCZ"]),
    tuple("AD-LON","${params.clump_base}/AD_LON/AD-LON",["AD","LON"]),
    tuple("SCZ-LON","${params.clump_base}/SCZ_LON/SCZ-LON",["SCZ","LON"]),
    tuple("AD-SCZ-LON","${params.clump_base}/AD-SCZ-LON/AD-SCZ",["AD","SCZ","LON"])
  ) | PREP_4_FUMA

  Channel.of(
    tuple("AD","${params.clump_base}/AD_SCZ/AD-SCZ"),
    tuple("SCZ","${params.clump_base}/AD_SCZ/AD-SCZ"),
    tuple("AD","${params.clump_base}/AD_LON/AD-LON"),
    tuple("LON","${params.clump_base}/AD_LON/AD-LON"),
    tuple("SCZ","${params.clump_base}/SCZ_LON/SCZ-LON"),
    tuple("LON","${params.clump_base}/SCZ_LON/SCZ-LON"),
    tuple("AD","${params.clump_base}/AD-SCZ-LON/AD-SCZ"),
    tuple("SCZ","${params.clump_base}/AD-SCZ-LON/AD-SCZ"),
    tuple("LON","${params.clump_base}/AD-SCZ-LON/AD-SCZ")
  ) | FUMA_POST_DIRS
}
