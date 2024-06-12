import scanpy as sc
import pandas as pd
import celltypist
from celltypist import models

#Here we can also source a script where we set our directory paths. Please change
#the file accordingly
with open('functions/python_directories.py') as f:
    exec(f.read())

models.download_models(force_update = True, model=['Immune_All_Low.pkl', 'Immune_All_High.pkl'])
model = models.Model.load()

m_adata = sc.read_h5ad(f'{vDir}/data/0_m_anndata.h5ad')
t_adata = sc.read_h5ad(f'{vDir}/data/0_t_anndata.h5ad')
b_adata = sc.read_h5ad(f'{vDir}/data/0_b_anndata.h5ad')

def predict_ct(data, prefix, model = model, majority_voting = False, mode = 'best match'):
    data.obsm['X_umap'] = data.obsm['X_UMAP']
    data.obsm['X_pca'] = data.obsm['X_PCA']
    data.var_names = data.var['human_symbol']
    predictions = celltypist.annotate(data, model = model, majority_voting = majority_voting, mode = mode)
    pdata = predictions.to_adata()
    out = pdata.obs
    out.to_csv(f'{plotsDir}{prefix}_res_celltypist.csv')

dic = {'myeloid': m_adata, 'lymphoid': t_adata, 'b': b_adata}

[predict_ct(data = value, prefix = key) for (key,value) in dic.items()]

