import scanpy as sc
import sys

adata = sc.read_h5ad('brain_data.h5ad')
print(f'Shape: {adata.n_obs} cells x {adata.n_vars} genes')
print(f'obs columns: {list(adata.obs.columns)}')
for col in ['cell_type', 'celltype', 'tissue', 'tissue_general', 'disease', 'organism']:
    if col in adata.obs.columns:
        vc = adata.obs[col].value_counts()
        print(f'\n{col} value counts (top 20):')
        print(vc.head(20).to_string())
