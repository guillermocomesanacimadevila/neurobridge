#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.pybin = "python3"
params.smr_alpha = 0.05
params.heidi_alpha = 0.01
params.gene_list = ""
params.pheno_pair = "AD_SCZ"
params.pheno1 = "AD"
params.pheno2 = "SCZ"
params.smr_base = "${workflow.launchDir}/outputs/SMR/res"
params.out_base = "${workflow.launchDir}/outputs/SMR/res"
params.matrix_py = "${workflow.launchDir}/src/sMR/smr_fdr_matrix.py"
params.hits_py = "${workflow.launchDir}/src/sMR/hits.py"

// need a way to include locus coordinates (maybe from conjFDR stuff)
// input: gene lists -> each locus (1 gene list)
// out_dir (after the matrix script) -> outputs/SMR/${pheno_pair}/${pheno_id}/${gene_list[0]}_{gene_list[1]}_{gene_list[2]}_..._.tsv
// out_dir (after hits.py) -> outputs/SMR/${pheno_pair}/overlap/${gene_list[0]}_{gene_list[1]}_{gene_list[2]}_..._.tsv
// pheno1_id
// pheno2_id
// gene list
// p_heidi cutoff
// p_smr cutt off

// nextflow run nf/smr_overlap.nf \
//   --pheno_pair AD_SCZ \
//   --pheno1 AD \
//   --pheno2 SCZ \
//   --gene_list CLU,PBK,EPHX2,CHRNA2,PNOC,CCDC25

process FDR_AND_FILTER {
  tag "${pheno_id}"

  publishDir "${params.out_base}/${pheno_pair}/${pheno_id}", mode:'copy', saveAs: { fn ->
    fn.endsWith(".tsv") ? "${geneset}.tsv" : fn
  }

  input:
  tuple val(pheno_pair), val(pheno_id), val(geneset), val(gene_list), path(infile)

  output:
  tuple val(pheno_pair), val(geneset), val(pheno_id), path("${pheno_id}.${geneset}.tsv")

  script:
  """
  set -euo pipefail
  ${params.pybin} ${params.matrix_py} \\
    --in ${infile} \\
    --out_dir . \\
    --pheno_id ${geneset} \\
    --genes "${gene_list}" \\
    --pcol p_SMR_multi \\
    --alpha ${params.smr_alpha} \\
    --pheidi_col p_HEIDI \\
    --pheidi_min ${params.heidi_alpha}

  mv ${geneset}.tsv ${pheno_id}.${geneset}.tsv
  """
}

process OVERLAPPING_HITS {
  tag "${pheno_pair}"
  publishDir "${params.out_base}/${pheno_pair}/overlaps", mode:'copy'

  input:
  tuple val(pheno_pair), val(geneset), path(mat_a), path(mat_b)

  output:
  path("${geneset}.tsv")

  script:
  """
  set -euo pipefail
  ${params.pybin} ${params.hits_py} \\
    --mat_a ${mat_a} \\
    --mat_b ${mat_b} \\
    --out ${geneset}.tsv \\
    --trait_a ${params.pheno1} \\
    --trait_b ${params.pheno2} \\
    --alpha ${params.smr_alpha}
  """
}

workflow {
  if( !params.gene_list?.trim() ) error "need --gene_list \"GENE1,GENE2,...\""
  def genes = params.gene_list.split(',').collect{ it.trim() }.findAll{ it }
  def geneset = genes.join('_')
  def in1 = file("${params.smr_base}/${params.pheno1}/trait_eSMR.merged.tsv")
  def in2 = file("${params.smr_base}/${params.pheno2}/trait_eSMR.merged.tsv")
  Channel.of(
    tuple(params.pheno_pair, params.pheno1, geneset, params.gene_list, in1),
    tuple(params.pheno_pair, params.pheno2, geneset, params.gene_list, in2)
  ) | FDR_AND_FILTER | set { mats }
  def a = mats.filter{ it[2] == params.pheno1 }.map{ it[3] }
  def b = mats.filter{ it[2] == params.pheno2 }.map{ it[3] }

  a.combine(b).map{ ma, mb -> tuple(params.pheno_pair, geneset, ma, mb) } | OVERLAPPING_HITS
}
