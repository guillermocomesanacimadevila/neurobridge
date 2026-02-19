import pandas as pd
from statsmodels.stats.multitest import multipletests

pd.set_option("display.float_format", "{:.10f}".format)
genes = [
    "ACP2","MADD","ADAM10","MINDY2","CLU","EPHX2","PBK","CHRNA2","PNOC","CCDC25",
    "RNF111","ALDH1A2","AQP9","CCNB2","GCNT3","LDHAL6B","LIPC","SLTM","SPI1","AGBL2",
    "ARFGAP2","ARHGAP1","ATG13","C11orf49","C1QTNF4","CELF1","CKAP5","DDB2","F2",
    "FAM180B","FNBP4","LRP4","MDK","MTCH2","MYBPC3","NDUFS3","NR1H3","PACSIN3",
    "PSMC3","PTPMT1","PTPRJ","RAPSN","SLC39A13","ZNF408","NUP160"
]

ad = pd.read_csv("ad_smr.tsv", sep="\t")
scz = pd.read_csv("scz_smr.tsv", sep="\t")

def add_q_gene_fdr(df, qtl):
    sub = df[df["qtl_name"] == qtl].copy()
    sub["q_SMR_multi_gene_fdr"] = multipletests(sub["p_SMR_multi"].astype(float), method="fdr_bh")[1]
    return sub

for qtl in ["eQTL_eQTLGen", "eQTL_BrainMeta"]:
    ad_full = add_q_gene_fdr(ad, qtl)
    scz_full = add_q_gene_fdr(scz, qtl)

    ad_full.to_csv(f"AD_{qtl}_allgenes.tsv", sep="\t", index=False)
    scz_full.to_csv(f"SCZ_{qtl}_allgenes.tsv", sep="\t", index=False)

print(len(ad["qtl_name"].unique()))
print(len(scz["qtl_name"].unique()))
for qtl in ad["qtl_name"].unique():
    ad_sub = ad[ad["qtl_name"] == qtl].copy()
    scz_sub = scz[scz["qtl_name"] == qtl].copy()
    ad_sub["q_SMR_multi"] = multipletests(ad_sub["p_SMR_multi"].astype(float), method="fdr_bh")[1]
    scz_sub["q_SMR_multi"] = multipletests(scz_sub["p_SMR_multi"].astype(float), method="fdr_bh")[1]
    rows = []
    for gene in genes:
        ad_gene = ad_sub[ad_sub["index"].str.upper() == gene][["index","p_SMR_multi","q_SMR_multi","p_HEIDI"]].copy()
        scz_gene = scz_sub[scz_sub["index"].str.upper() == gene][["index","p_SMR_multi","q_SMR_multi","p_HEIDI"]].copy()
        ad_gene["trait"] = "AD"
        scz_gene["trait"] = "SCZ"
        rows.append(ad_gene)
        rows.append(scz_gene)
    out = pd.concat(rows, ignore_index=True)
    def keep_gene(df):
        if df.shape[0] != 2:
            return False
        cond1 = (df["p_HEIDI"].astype(float) > 0).any() # 0.01
        cond2 = (df["q_SMR_multi"].astype(float) < 0.05).all()
        return cond1 and cond2
    out = (
        out.groupby(out["index"].str.upper(), group_keys=False)
           .filter(keep_gene)
    )
    print(f"\n=== {qtl} ===")
    print(out.sort_values(["index","trait"]))
