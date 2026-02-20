import cellxgene_census

# make another scripts that combined all gene IDs from sMR outputs
# CLI args ->
# combined gene file
# cell types which interested
# tissues interested
# normal OR disease -> if disease -> which one
# out_dir
# pheno1_prefix
# pheno2_prefix

gene = "SPI1"
cell_type = "microglial cell"
census = cellxgene_census.open_soma(census_version="2025-11-08")
obs_filter = f'tissue_general == "brain" and cell_type == "{cell_type}" and disease == "normal"'
var_filter = f'feature_name == "{gene}"'
adata = cellxgene_census.get_anndata(
    census,
    organism="homo_sapiens",
    obs_value_filter=obs_filter,
    var_value_filter=var_filter
)

print(adata)
census.close()
