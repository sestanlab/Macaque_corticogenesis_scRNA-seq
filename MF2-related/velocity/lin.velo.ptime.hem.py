'''import required modules'''
import scvelo as scv
import pandas as pd
import numpy as np
import seaborn as sns

lin = "Hem"

out_h5ad = "./load_files/velo_" + lin + "_org_addptime.h5ad"
adata = scv.read("./load_files/velo_" + lin + "_org.h5ad")



scv.tl.latent_time(adata)
scv.settings.set_figure_params(dpi_save = 300, transparent = True)
scv.pl.scatter(adata, save = "scVelo_" + lin + "_latent_time.png", color='latent_time', color_map='gnuplot', size=30, colorbar=True, figsize = (7,7))

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
sns.set(font_scale=0.5)
scv.settings.set_figure_params(dpi_save = 300, transparent = False)

scv.pl.heatmap(adata, save = "scVelo_" + lin + "_gene_cascade_full.png",var_names=top_genes, tkey='latent_time', n_convolve=100, col_color='cluster', figsize = (8, 40), yticklabels=True)



## Store the dataset
adata.write(out_h5ad)



##------------------------------------------------------
## Get ordered and fitted data [Same order in the heatmap]
xkey = 'Ms'
tkey = 'latent_time'
n_convolve = 100

var_names = top_genes[:300]
add_genes = ["NKX2-1", "GNRH1", "DLX1", "GAD1", "LHX6"]
for i in add_genes:
    if i not in var_names and i in adata.var_names:
            var_names = var_names.append(pd.Index([i]))


time = adata.obs[tkey].values
time = time[np.isfinite(time)]

df = pd.DataFrame(adata[:, var_names].layers[xkey][np.argsort(time)], columns=var_names, index = adata.obs.index[np.argsort(time)])

if n_convolve is not None:
    weights = np.ones(n_convolve) / n_convolve
    for i, gene in enumerate(var_names):
        df[gene] = np.convolve(df[gene].values, weights, mode='same')

max_sort = np.argsort(np.argmax(df.values, axis=0))
df = pd.DataFrame(df.values[:, max_sort], columns=df.columns[max_sort], index = df.index)
df.to_csv(path_or_buf = "./load_files/Heatdata." + lin + ".txt", sep = "\t")



## Write out latent time
time = adata.obs[['latent_time', 'cluster']]
time.to_csv(path_or_buf = "./load_files/Latentime." + lin + ".txt", sep = "\t")


