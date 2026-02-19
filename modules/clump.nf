#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.pybin = "python3"
params.ref_bfile = "${workflow.launchDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL"
params.clump_out = "outputs/clumping"

process CLUMP_PAIRWISE {
  tag "${pair_id}"

  input:
  tuple val(pair_id),
        path(hits),
        val(gwas1),
        val(gwas2),
        val(p1),
        val(p2),
        val(out_dir)

  output:
  path("${pair_id}_clumping_prep.done")

  script:
  """
  set -e
  mkdir -p ${out_dir}
  ${params.pybin} ${workflow.launchDir}/src/ld-clump/clump.py \\
    --mode pairwise \\
    --hits ${hits} \\
    --out-dir ${out_dir} \\
    --pheno1 ${p1} \\
    --pheno2 ${p2} \\
    --ref-bfile ${params.ref_bfile} \\
    --score-col conj_fdr \\
    --clump-r2 0.6 \\
    --clump-kb 1000 \\
    --merge-kb 500 \\
    --ld-r2 0.1 \\
    --ld-kb 250 \\
    --sumstats ${gwas1} ${gwas2} \\
    --sumstats-names ${p1} ${p2}

  touch ${pair_id}_clumping_prep.done
  """
}

process CLUMP_TRIPLE_OVERLAP {
  tag "${triple_id}"

  input:
  tuple val(triple_id),
        path(hits),
        val(gwas1),
        val(gwas2),
        val(gwas3),
        val(p1),
        val(p2),
        val(p3),
        val(triple_col1),
        val(triple_col2),
        val(out_dir)

  output:
  path("${triple_id}_clumping_prep.done")

  script:
  """
  set -e
  mkdir -p ${out_dir}
  ${params.pybin} ${workflow.launchDir}/src/ld-clump/clump.py \\
    --mode triple \\
    --hits ${hits} \\
    --out-dir ${out_dir} \\
    --pheno1 ${p1} \\
    --pheno2 ${p2} \\
    --ref-bfile ${params.ref_bfile} \\
    --clump-r2 0.6 \\
    --clump-kb 1000 \\
    --merge-kb 500 \\
    --ld-r2 0.1 \\
    --ld-kb 250 \\
    --triple-col1 ${triple_col1} \\
    --triple-col2 ${triple_col2} \\
    --triple-how max \\
    --sumstats ${gwas1} ${gwas2} ${gwas3} \\
    --sumstats-names ${p1} ${p2} ${p3}
    
  touch ${triple_id}_clumping_prep.done
  """
}

workflow {
  def base = "${workflow.launchDir}/data"
  def conj = "${workflow.launchDir}/outputs/conjFDR"
  def outb = "${workflow.launchDir}/${params.clump_out}"
  def ad  = file("${base}/AD/post-qc/AD.ldsc_ready_neff.tsv").toAbsolutePath().toString()
  def scz = file("${base}/SCZ/post-qc/SCZ.ldsc_ready_neff.tsv").toAbsolutePath().toString()
  def lon = file("${base}/LON/post-qc/LON.ldsc_ready_neff.tsv").toAbsolutePath().toString()

  pairwise_info = Channel.of(
    tuple(
      "AD_SCZ",
      file("${conj}/AD_SCZ/AD_SCZ_shared_hits.tsv"),
      ad,
      scz,
      "AD",
      "SCZ",
      "${outb}/AD_SCZ"
    ),
    tuple(
      "AD_LON",
      file("${conj}/AD_LON/AD_LON_shared_hits.tsv"),
      ad,
      lon,
      "AD",
      "LON",
      "${outb}/AD_LON"
    ),
    tuple(
      "SCZ_LON",
      file("${conj}/SCZ_LON/SCZ_LON_shared_hits.tsv"),
      scz,
      lon,
      "SCZ",
      "LON",
      "${outb}/SCZ_LON"
    )
  )
  triple_info = Channel.of(
    tuple(
      "AD-SCZ-LON",
      file("${conj}/AD-SCZ-LON/overlap_merged.tsv"),
      ad,
      scz,
      lon,
      "AD",
      "SCZ",
      "LON",
      "conj_fdr_AD_LON",
      "conj_fdr_SCZ_LON",
      "${outb}/AD-SCZ-LON"
    )
  )
  CLUMP_PAIRWISE(pairwise_info)
  CLUMP_TRIPLE_OVERLAP(triple_info)
}
