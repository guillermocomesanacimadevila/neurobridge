# nf-core/neurobridge (xQtl Colocalization): Usage

```bash
python smr_hits.py \
  --pheno1_id AD \
  --pheno2_id SCZ \
  --pheno1_df ../../dat/AD/AD_trait_sSMR.merged.tsv \
  --pheno2_df ../../dat/SCZ/SCZ_trait_sSMR.merged.tsv \
  --xqtl sQTL \
  --chr 15 \
  --start 58534174 \
  --end 59534174 \
  --outdir results
```

```bash
python ld_smr_susie.py \
  --pheno1_id AD \
  --pheno2_id SCZ \
  --xqtl sQTL \
  --res_dir results \
  --chr 15 \
  --start 58534174 \
  --end 59534174 \
  --susie_dir /Users/c24102394/Desktop/neurobridge/outputs/susie/res \
  --ref_bfile /Users/c24102394/Desktop/neurobridge/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL \
  --outdir results/ld_results
```

```bash
python preprocess_xQTL_sumstats_for_coloc.py \
  outputs/xqtl_coloc_ready \
  sQTL \
  2865 \
  BrainMeta
```

```bash
python prep_gwas_and_xqtl_for_coloc.py \
  --pheno1_id AD \
  --pheno2_id SCZ \
  --ld_dir results/ld_results \
  --qtl_type sQTL \
  --qtl_dataset BrainMeta \
  --gwas_loci_dir /Users/c24102394/Desktop/neurobridge/outputs/defined_loci \
  --qtl_loci_dir outputs/xqtl_coloc_ready \
  --chr 15 \
  --smr_dir results \
  --start_locus 58534174 \
  --end_locus 59534174 \
  --window 250000 \
  --out_dir results/coloc_inputs
```

