#!/usr/bin/env python3
# https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_quick_start.html
import numpy as np
import pandas as pd
import cellxgene_census
import scanpy as sc
# scRNA env

GENE_ID = "ENSG00000066336"
TISSUE_GENERAL = "brain"
with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",
        measurement_name="RNA",
        var_value_filter=f"feature_id == '{GENE_ID}'",
        obs_value_filter=f"tissue_general == '{TISSUE_GENERAL}'",
        column_names={
            "obs": ["cell_type", "tissue", "tissue_general", "disease"],
            "var": ["feature_id", "feature_name"]
        },
    )

print(adata)
print("Gene:", adata.var["feature_name"].values[0])
