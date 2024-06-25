import numpy as np
import pandas as pd
import scanpy as sc

import plotnine as p9

import liana as li
import decoupler as dc
import omnipath as op
import corneto as cn
import gurobipy

from plotnine import *

adata = sc.read_h5ad(filename=f'{dataDir}/3_full_anndata.h5ad')

dge_res=pd.read_csv(f'{resDir}/5_dge_full_dataset.csv')
dge_res['coef'] = pd.Categorical(dge_res.coef)

comparisons = dge_res['coef'].unique()[::-1]

def differential_interactions(comp, org, exp):
    dge = dge_res.query('coef == @comp & organ == @org & experiment == @exp')
    dge = dge.set_index('rn')
    dge = dge.rename(columns={'t':'t_stat'})
    cond = comp[0:comp.find('_vs_')]
    subset = adata[(adata.obs['treatment.agg'] == cond) & (adata.obs['organ'] == org) & (adata.obs['experiment'] == exp)].copy()
    df = li.multi.df_to_lr(subset, 
                           dea_df=dge, 
                           groupby='celltype', 
                           resource_name='mouseconsensus',
                           stat_keys=['t_stat', 'P.Value', 'adj.P.Val', 'B'],
                           complex_col='t_stat',
                           expr_prop=0.1,
                           use_raw=False, 
                           verbose=True,
                           return_all_lrs=False)
    return(df)

organ = adata.obs['organ'].copy().unique()
experiment = adata.obs['experiment'].copy().unique()

lr_res={f'{i}-{j}-{m}':differential_interactions(i, j, m) for i in comparisons for j in organ for m in experiment}

lr_res = pd.concat(lr_res.values(), keys = lr_res.keys())
lr_res = lr_res.reset_index(level=0)
lr_res = lr_res.rename({'level_0':'analysis'}, axis=1)
lr_res[['comparison', 'organ', 'experiment']] = lr_res['analysis'].str.split('-', expand = True)
lr_res[['contrast', 'intercept']] = lr_res['comparison'].str.split('_vs_', expand = True)

lr_res.to_csv(f'{resDir}/9_liana_differntial_interactions.csv', index=False)
lr_res = pd.read_csv(f'{resDir}/9_liana_differntial_interactions.csv')

lr_res = lr_res.sort_values(['interaction_t_stat', 'analysis'], ascending=[False,True])

lr_res.hist(bins=50, by = 'analysis', column = 'interaction_t_stat')

for i in lr_res['analysis'].unique():
    df = lr_res.query('analysis == @i')
    exp = df['experiment'].unique()[0]
    org = df['organ'].unique()[0]
    con = df['contrast'].unique()[0]
    inter = df['intercept'].unique()[0]

    f = (li.pl.tileplot(liana_res=df,
                       fill='expr',
                       label='adj.P.Val',
                       label_fun = lambda x: '*' if x < 0.05 else np.nan,
                       top_n=15,
                       orderby = 'interaction_t_stat',
                       orderby_ascending = False,
                       orderby_absolute = False,
                       source_title='Ligand',
                       target_title='Receptor',
                       figure_size=(11,10)
                      ) +
         p9.ggtitle(f'{con} vs {inter} | {org} | {exp}'))
    f.save(f'{plotsDir}/{con}_vs_{inter}_{org}_{exp}.pdf')

#Set the comparison as one of the following:
lr_res['comparison'].unique()
#So.. HDAC_cKO_vs_NoT, HDAC_cKO_vs_HDAC_WT or HDAC_WT_vs_NoT
#depending on your interest

comp='HDAC_cKO_vs_HDAC_WT'

#Now which experiment group? HDAC1 or HDAC2
exp='HDAC2'

t = lr_res.query('comparison == @comp & experiment == @exp')

t_n_inter=t.query('`interaction_adj.P.Val`<0.05').groupby(['source', 'target', 'organ'], as_index=False, observed=True)['interaction_props'].agg(count='count')

p=(ggplot(t_n_inter)+
   geom_tile(aes('source','target', fill='count'))+
   facet_grid('~organ', scales = 'free', space = 'free')+
   scale_fill_gradient(low = '#fee0d2', high = '#de2d26')+
   theme_bw()+
   theme(axis_text_x=element_text(rotation = 90)))

p.save(f'{plotsDir}/9_n_dge_interactions_HDAC_cKO.pdf', width=10, height = 10)

ligand_producing_cells=['T_CD45.1_2', 'T_gd_3', 'Keratinocytes_4', 'Fibroblasts_2']
receptor_producing_cells=['T_CD45.1_3', 'Fibroblasts_2', 'T_gd_3', 'Keratinocytes_3']
cut_plot = (li.pl.dotplot(liana_res=t,
                     colour='interaction_t_stat',
                     size='ligand_adj.P.Val',
                     inverse_size=True,
                     source_labels=ligand_producing_cells,
                     target_labels=receptor_producing_cells,
                     orderby='interaction_t_stat',
                     orderby_ascending=False,
                     orderby_absolute=True,
                     top_n=20,
                     size_range=(0.5, 4)
                     )

# customize plot

    + p9.theme_bw(base_size=14)
    # fill cmap blue to red, with 0 the middle
    + p9.scale_color_gradient2(low = 'blue', high = 'red', mid = 'white')
    # rotate x
    + p9.theme(axis_text_x=p9.element_text(angle=90), figure_size=(11, 6))

           )
           
cut_plot.save(f'{plotsDir}/9_top_differential_interactions_{comp}_{exp}.pdf')

def select_top_n(d, n=None):
    d = dict(sorted(d.items(), key=lambda item: abs(item[1]), reverse=True))
    return {k: v for i, (k, v) in enumerate(d.items()) if i < n}    

source_label = 'T_CD45.1_2'
target_label = 'Keratinocytes_3'

org='Skin'

# NOTE: We sort by the absolute value of the interaction stat
lr_stats=lr_res.query('`interaction_adj.P.Val` < 0.05 & comparison == @comp & experiment == @exp & organ == @org & source == @source_label & target == @target_label').copy()

lr_stats = lr_stats.sort_values('interaction_t_stat', ascending=False, key=abs)

lr_dict = lr_stats.set_index('receptor')['interaction_t_stat'].to_dict()
input_scores = select_top_n(lr_dict, n=15)

dge_res['treatment']

dge=dge_res.query('organ == @org & experiment == @exp &  treatment=="KO_vs_WT"')

# First, let's transform the DEA statistics into a DF
# we will use these to estimate deregulated TF activity
dge_wide = dge[['celltype', 't', 'rn']].pivot(index='celltype', columns='rn', values='t')
dge_wide = dge_wide.fillna(0)

net = dc.get_collectri(organism='mouse')
# Run Enrichment Analysis
estimates, pvals = dc.run_ulm(mat=dge_wide, net=net)
estimates.T.sort_values(target_label, key=abs, ascending=False).head()


tf_data = estimates.copy() 
tf_dict = tf_data.loc[target_label]
tf_dict = tf_dict[pvals.loc[target_label]<0.05].to_dict()
output_scores = select_top_n(tf_dict, n=10)


ppis = op.interactions.OmniPath().get(genesymbols = True, organism = 'mouse')

ppis['mor'] = ppis['is_stimulation'].astype(int) - ppis['is_inhibition'].astype(int)
ppis = ppis[(ppis['mor'] != 0) & (ppis['curation_effort'] >= 5) & ppis['consensus_direction']]

input_pkn = ppis[['source_genesymbol', 'mor', 'target_genesymbol']]
input_pkn.columns = ['source', 'mor', 'target']


prior_graph = li.mt.build_prior_network(input_pkn, input_scores, output_scores, verbose=True)

temp = adata[(adata.obs['celltype'] == target_label) & (adata.obs['experiment'] == 'HDAC2') & (adata.obs['organ'] == 'Skin') & (adata.obs['treatment.agg'] == 'HDAC_WT')].copy()

node_weights = pd.DataFrame(temp.X.getnnz(axis=0) / temp.n_obs, index=temp.var_names)
node_weights = node_weights.rename(columns={0: 'props'})
node_weights = node_weights['props'].to_dict()


df_res, problem = li.mt.find_causalnet(
    prior_graph,
    input_scores,
    output_scores,
    node_weights,
    # penalize (max_penalty) nodes with counts in less than 0.1 of the cells
    node_cutoff=0.1,
    max_penalty=1,
    # the penaly of those in > 0.1 prop of cells set to:
    min_penalty=0.01,
    edge_penalty=0.1,
    verbose=True,
    solver='scipy' # 'scipy' is available by default, but results in suboptimal solutions
    )
    
p = cn.methods.carnival.visualize_network(df_res)
p
p.render(filename = f'{plotsDir}/diff_network_KC_TGD_HDAC_cKO_vs_WT')
