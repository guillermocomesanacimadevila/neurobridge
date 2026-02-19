#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path
import argparse

# steps
# load res from FUMA cell-type
# load respective magma_celltype_step1.txt
# FDR afjusted p col = P.adj.pds
# for each given overlapping dataset -> check overlapping cell type -> if FDR < 0.05 on both
# Save row and rename cols adjusted to specific trait to outputs/sc-Enrichment/res.tsv

def get_cell_type_overlap(trait1: str,
                          trait2: str,
                          out_dir: str,
                          fdr_thresh: float,
                          trait1_prefix: str,
                          trait2_prefix: str):
    mode = "intersection" # union
    trait1 = pd.read_csv(trait1, sep="\t")
    trait2 = pd.read_csv(trait2, sep="\t")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ds1 = trait1["Dataset"].unique().tolist()
    ds2 = trait2["Dataset"].unique().tolist()
    shared = []
    for i in ds1:
        for j in ds2:
            if i == j:
                shared.append(i)
    t1 = trait1[trait1["Dataset"].isin(shared)]
    t2 = trait2[trait2["Dataset"].isin(shared)]
    t1["P.adj.pds"] = pd.to_numeric(t1["P.adj.pds"], errors="coerce")
    t2["P.adj.pds"] = pd.to_numeric(t2["P.adj.pds"], errors="coerce")
    t1 = t1[t1["P.adj.pds"] <= fdr_thresh]
    t2 = t2[t2["P.adj.pds"] <= fdr_thresh]
    keys1 = []
    keys2 = []
    for index, row in t1[["Dataset", "Cell_type"]].drop_duplicates().iterrows():
        keys1.append((row["Dataset"], row["Cell_type"]))
    for index, row in t2[["Dataset", "Cell_type"]].drop_duplicates().iterrows():
        keys2.append((row["Dataset"], row["Cell_type"]))
    selected_keys = (
        [k for k in keys1 if k in keys2]
        if mode == "intersection"
        else list(set(keys1 + keys2))
    )
    t1["__key__"] = list(zip(t1["Dataset"], t1["Cell_type"]))
    t1 = t1[t1["__key__"].isin(selected_keys)].drop(columns="__key__")
    t2["__key__"] = list(zip(t2["Dataset"], t2["Cell_type"]))
    t2 = t2[t2["__key__"].isin(selected_keys)].drop(columns="__key__")
    df = pd.merge(
        t1,
        t2,
        on=["Dataset", "Cell_type"],
        how="inner" if mode == "intersection" else "outer",
        suffixes=(trait1_prefix, trait2_prefix)
    )
    df[f"sig{trait1_prefix}"] = df.get(f"P.adj.pds{trait1_prefix}").notna()
    df[f"sig{trait2_prefix}"] = df.get(f"P.adj.pds{trait2_prefix}").notna()
    out_file = out_dir / f"{trait1_prefix}_{trait2_prefix}_FUMA_sc_enrichment_{mode}.tsv"
    df.to_csv(out_file, sep="\t", index=False)
    return df["Cell_type"].unique()

# for each cell dataset
# for each cell type within each shared dataset
# if cell type present in df1 and df2 and padj < 0.05  keep

# /Users/c24102394/Desktop/FUMA/res/AD/magma_celltype_step1.txt
ad = "/Users/c24102394/Desktop/FUMA/res/AD/magma_celltype_step1.txt"
scz = "/Users/c24102394/Desktop/FUMA/res/SCZ/magma_celltype_step1.txt"
print(get_cell_type_overlap(ad, scz, "/Users/c24102394/Desktop/FUMA/res", 0.05, "AD", "SCZ"))