'''import required modules'''
import scvelo as scv
import pandas as pd
import numpy as np
import seaborn as sns

lin = "Hem"

out_h5ad = "./load_files/velo_" + lin + "_org.h5ad"
adata = scv.read("./load_files/velo." + lin + ".h5ad")


## Add colors
cls_use = list(adata.obs['cluster'].astype('category').cat.categories)
col_anno = pd.read_csv("./load_files/cls.col.txt", sep = "\t", header = None, names=["cluster", "color"])
col_use = []
for i in cls_use:
    if i in col_anno["cluster"].tolist():
        new_col = col_anno["color"].tolist()[col_anno["cluster"].tolist().index(i)]
        col_use.append(new_col)
adata.uns['cluster_colors']=col_use
print(adata.uns['cluster_colors'])



scv.pp.filter_and_normalize(adata, min_shared_counts=8, n_top_genes=1000)
scv.pp.moments(adata, n_pcs=25, n_neighbors=20)
scv.tl.recover_dynamics(adata, max_iter=400, return_model=True, n_jobs = 12)
scv.tl.velocity(adata, mode='dynamical')
##scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", save = "scVelo_" + lin + "_stream.png", dpi=300, color = "cluster", size = 17, fontsize = 0, title = "", alpha = 1, legend_loc='none', linewidth = 0.5, density = 0.8, figsize = (4,4))
scv.pl.velocity_embedding(adata, basis='umap', save = "scVelo_" + lin + "_arrow.png", dpi=300, arrow_length=2, arrow_size=1.5, color = "cluster", fontsize = 0, title = "", alpha = 1, legend_loc='none', size = 12, figsize = (4,4))
scv.pl.velocity_embedding_grid(adata, basis='umap',save = "scVelo_" + lin + "_grid.png", dpi=300,  arrow_length = 2, arrow_size = 1.5, color = "cluster", fontsize = 0, title = "", alpha = 1, legend_loc='none', size = 17, density = 0.8, figsize = (4,4))





##------------------------------------------------------
## Store the dataset
## There is a renaming error: the problem is described here
## https://github.com/theislab/scvelo/issues/255
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata.write(out_h5ad)

