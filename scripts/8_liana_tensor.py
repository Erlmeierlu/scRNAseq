import pandas as pd
import scanpy as sc
import plotnine as p9

import liana as li
import cell2cell as c2c
import decoupler as dc # needed for pathway enrichment

import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict

#Directories..
with open('functions/python_directories.py') as f:
    exec(f.read())
    
adata = sc.read_h5ad(filename=f'{dataDir}/7_full_anndata.h5ad')

#This script follows a liana+ tutorial. Most things are explained there:
#https://liana-py.readthedocs.io/en/latest/notebooks/liana_c2c.html

#For sample-wise analysis we pass the sample_treat column as the sample 
#key. The sample treat column is a combination of the sample, and the assigned
#treatment group. (e.g. fLN_40B2___HDAC_WT)
sample_key = 'sample_treat.agg'
condition_key = 'treatment.agg'
groupby = 'celltype'

#here we need to remove some sample_treatment combination with less than 3 occurences
#Otherwise the analysis would fail..
samples_used = adata.obs[adata.obs[sample_key].map(adata.obs[sample_key].value_counts()) > 2][sample_key].cat.remove_unused_categories().unique()

adata = adata[adata.obs[sample_key].isin(samples_used)].copy()

#this will run for a long time
li.mt.rank_aggregate.by_sample(
    adata,
    groupby=groupby,
    min_cells=0,
    sample_key=sample_key, # sample key by which we which to loop
    use_raw=False,
    resource_name='mouseconsensus',
    verbose=True, # use 'full' to show all verbose information
    n_perms=200, # reduce permutations for speed
    return_all_lrs=True, # return all LR values
    )
    
adata.uns["liana_res"] = pd.concat(adata.uns["liana_res"].values(), keys = adata.uns["liana_res"].keys())
adata.uns["liana_res"] = adata.uns["liana_res"].reset_index(level=0).rename({'level_0':'sample_treat'}, axis=1)

tensor = li.multi.to_tensor_c2c(adata,
                                sample_key=sample_key,
                                score_key='magnitude_rank', # can be any score from liana
                                how='outer_cells' # how to join the samples
                                )

#Exporting Tensor & Loading it
c2c.io.export_variable_with_pickle(tensor, f'{dataDir}/8_tensor_liana.pkl')
tensor=c2c.io.load_variable_with_pickle(f'{dataDir}/8_tensor_liana.pkl')

context_dict = adata.obs[[sample_key, condition_key]].drop_duplicates()
context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
context_dict = defaultdict(lambda: 'Unknown', context_dict)

tensor_meta = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                  metadata_dicts=[context_dict, None, None, None],
                                                  fill_with_order_elements=True
                                                  )

#This will also run for some time!
tensor = c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                                    tensor_meta,
                                                    copy_tensor=True, # Whether to output a new tensor or modifying the original
                                                    rank=None, # Number of factors to perform the factorization. If None, it is automatically determined by an elbow analysis. Here, it was precomuputed.
                                                    tf_optimization='regular', # To define how robust we want the analysis to be.
                                                    random_state=None, # Random seed for reproducibility
                                                    device='cpu', # Device to use. If using GPU and PyTorch, use 'cuda'. For CPU use 'cpu'
                                                    elbow_metric='error', # Metric to use in the elbow analysis.
                                                    smooth_elbow=False, # Whether smoothing the metric of the elbow analysis.
                                                    upper_rank=25, # Max number of factors to try in the elbow analysis
                                                    tf_init='random', # Initialization method of the tensor factorization
                                                    tf_svd='numpy_svd', # Type of SVD to use if the initialization is 'svd'
                                                    cmaps=None, # Color palettes to use in color each of the dimensions. Must be a list of palettes.
                                                    sample_col='Element', # Columns containing the elements in the tensor metadata
                                                    group_col='Category', # Columns containing the major groups in the tensor metadata
                                                    output_fig=False, # Whether to output the figures. If False, figures won't be saved a files if a folder was passed in output_folder.
                                                    )
                                                    
factors, axes = c2c.plotting.tensor_factors_plot(interaction_tensor=tensor,
                                                 metadata = tensor_meta, # This is the metadata for each dimension
                                                 sample_col='Element',
                                                 group_col='Category',
                                                 meta_cmaps = ['viridis', 'Dark2_r', 'tab20', 'tab20'],
                                                 fontsize=10, # Font size of the figures generated
                                                 )

factors = tensor.factors
factors.keys()

#From this plot we can determine which factors we want to look at in detail. 
#look at the liana+ tutorial what they did..
_ = c2c.plotting.context_boxplot(context_loadings=factors['Contexts'],
                                 metadict=context_dict,
                                 nrows=2,
                                 figsize=(8, 6),
                                 statistical_test='t-test_ind',
                                 pval_correction='fdr_bh',
                                 cmap='plasma',
                                 verbose=False,
                                )

#Here, also factor 6 was imortant.. Closer look..            
c2c.plotting.ccc_networks_plot(factors,
                               included_factors=['Factor 6'],
                               network_layout='circular',
                               ccc_threshold=0.05, # Only important communication
                               nrows=1,
                               panel_size=(8, 8), # This changes the size of each figure panel.
                              )
                              
lr_loadings = factors['Ligand-Receptor Pairs']
lr_loadings.sort_values("Factor 6", ascending=False).head(10)

# load PROGENy pathways
net = dc.get_progeny(organism='mouse', top=5000)

# load full list of ligand-receptor pairs
lr_pairs = li.resource.select_resource('mouseconsensus')
# generate ligand-receptor geneset
lr_progeny = li.rs.generate_lr_geneset(lr_pairs, net, lr_sep="^")
lr_progeny.head()

estimate, pvals = dc.run_ulm(lr_loadings.transpose(), lr_progeny, source="source", target="interaction", use_raw=False)
dc.plot_barplot(estimate, 'Factor 6', vertical=True, cmap='coolwarm', vmin=-7, vmax=7)
