import scanpy as sc
import liana as li
import numpy as np
import pandas as pd
import plotnine as p9

from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean, rank_aggregate
from plotnine import ggplot, geom_point, geom_tile, aes, facet_wrap, facet_grid, labs, theme_bw, theme, element_text, element_rect, scale_size_continuous, scale_fill_gradient, scale_color_cmap

#Directories..
with open('functions/python_directories.py') as f:
    exec(f.read())
    
adata = sc.read_h5ad(f'{dataDir}/3_full_anndata.h5ad')
adata.strings_to_categoricals()
adata.obs['treatment.agg'] = adata.obs['treatment.agg'].cat.reorder_categories(['NoT', 'HDAC_WT', 'HDAC_cKO', 'Undefined'])

levels = adata.obs['treatment.agg'].copy().sort_values().unique()
levels = levels[levels!='Undefined'].remove_unused_categories()
organ = adata.obs['organ'].copy().unique()
experiment = adata.obs['experiment'].copy().unique()

def analyse_conditions(obj, cond, org, exp):
    subset = obj[(obj.obs['treatment.agg'] == cond) & (obj.obs['organ'] == org) & (obj.obs['experiment'] == exp)].copy()
    df = rank_aggregate(subset, groupby='celltype', expr_prop=0.1, resource_name='mouseconsensus', use_raw=False, key_added = '', verbose = True, inplace = False)
    return(df)

res={f'{i}-{j}-{m}':analyse_conditions(adata, i, j, m) for i in levels for j in organ for m in experiment}

res = pd.concat(res.values(), keys = res.keys())
res = res.reset_index(level=0)

res = res.rename({'level_0':'analysis'}, axis=1)

res[['treatment', 'organ', 'experiment']] = res['analysis'].str.split('-', expand = True)

adata.uns['consensus_res'] = res

adata.uns['consensus_res']['organ'] = pd.Categorical(adata.uns['consensus_res'].organ)
adata.uns['consensus_res']['experiment'] = pd.Categorical(adata.uns['consensus_res'].experiment)
adata.uns['consensus_res']['treatment'] = pd.Categorical(adata.uns['consensus_res'].treatment)
adata.uns['consensus_res']['treatment'] = adata.uns['consensus_res']['treatment'].cat.reorder_categories(['NoT', 'HDAC_WT', 'HDAC_cKO'])


adata.obsm['X_umap'] = adata.obsm.pop('X_UMAP')
adata.obsm['X_pca'] = adata.obsm.pop('X_PCA')
adata.obsm['X_aligned'] = adata.obsm.pop('X_Aligned')

adata.write_h5ad(filename=f'{dataDir}/3_full_anndata.h5ad')
# adata = sc.read_h5ad(filename=f'{dataDir}/3_full_anndata.h5ad')

n_inter=adata.uns['consensus_res'].query('specificity_rank<0.05').groupby(['source', 'target', 'treatment', 'organ', 'experiment'], as_index=False, observed=True)['ligand_complex'].agg(count='count')

p=(ggplot(n_inter)+
   geom_tile(aes('source','target', fill='count'))+
   facet_grid('experiment+organ~treatment', scales = 'free', space = 'free')+
   scale_fill_gradient(low = '#fee0d2', high = '#de2d26')+
   theme_bw()+
   theme(axis_text_x=element_text(rotation = 90)))

p.save(f'{plotsDir}/consensus_n_interactions.pdf', width=15, height = 20)

a = adata.uns['consensus_res']
a['interaction'] = a['ligand_complex'] + ' -> ' + a['receptor_complex']
a = a.loc[a['cellphone_pvals'] <= 0.05]
a['-log10(p_vals)'] = -np.log10(a['cellphone_pvals'] + np.finfo(float).eps)


p = (ggplot(a.head(20), aes(x='target', y='interaction', colour='lr_means', size='-log10(p_vals)')) 
    + geom_point()
    + facet_grid('analysis~source')
     #+ scale_size_continuous(range=size_range)
     #+ scale_color_cmap(cmap)
     #+ labs(color=str.capitalize(colour),
     #       size=str.capitalize(size),
     #       y="Interactions (Ligand -> Receptor)",
     #       x="Target",
     #       title="Source")
     + theme_bw()
     + theme(legend_text=element_text(size=14),
             strip_background=element_rect(fill="white"),
             strip_text=element_text(size=15, colour="black"),
             axis_text_y=element_text(size=10, colour="black"),
             axis_title_y=element_text(colour="#808080", face="bold", size=15),
             axis_text_x=element_text(size=11, face="bold", angle=90),
             #figure_size=figure_size,
             plot_title=element_text(vjust=0, hjust=0.5, face="bold",
                                     colour="#808080", size=15)
             )
     )
p

