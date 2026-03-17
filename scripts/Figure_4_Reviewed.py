# %%
import anndata as ad
import pandas as pd
import numpy as np
import h5py
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import matplotlib.pyplot as plt

# %%
adata = ad.read_h5ad("data/wmb_sa2_anndata/sa2_logCPM_gmat.h5ad", backed='r')

# %%

meta = pd.read_csv("data/SI Table 2 Metadata table for all the 2.3 million nuclei in the snATAC-seq data.txt", 
                   sep="\t", index_col="CellID")

print(meta.shape)
print(meta['Subclass'].value_counts().head(20))

subclases_interes = ['L2/3 IT CTX', 'L5 IT CTX', 'L5 PT CTX', 'L6 IT CTX', 
                     'Pvalb', 'Sst', 'Vip', 'Lamp5']

for s in subclases_interes:
    matches = meta['Subclass'].str.contains(s, na=False).sum()
    print(f"{s}: {matches} células")

# %%
import h5py
import numpy as np
import pandas as pd


subclases_atac = list(subclases_map.values())
meta_filtrada = meta[meta['Subclass'].isin(subclases_atac)].copy()
inv_map = {v: k for k, v in subclases_map.items()}
meta_filtrada['subclass_short'] = meta_filtrada['Subclass'].map(inv_map)
genes_interes = ['Rpl22', 'Rpl22l1', 'Rpl7', 'Rpl7l1', 'Rps27', 'Rps27l']
f = h5py.File("data/sa2_logCPM_gmat.h5ad", 'r')


gene_names = pd.Index(f['var']['_index'][:].astype(str))
cell_names = pd.Index(f['obs']['_index'][:].astype(str))


cells_interes = meta_filtrada.index.intersection(cell_names)
print(f"Celulas de interes: {len(cells_interes)}")


gene_positions = np.array([gene_names.get_loc(g) for g in genes_interes])
cell_positions = np.where(cell_names.isin(cells_interes))[0]
print(f"Posiciones calculadas: {len(cell_positions)} celulas, {len(gene_positions)} genes")

# %%

indptr = f['X']['indptr'][:]
starts = indptr[cell_positions]
ends = indptr[cell_positions + 1]
total_entries = (ends - starts).sum()
print(f"Total entradas a leer: {total_entries:,}")
print(f"Memoria estimada: {total_entries * 8 / 1e9:.2f} GB")

# %%

all_starts = indptr[cell_positions]
all_ends = indptr[cell_positions + 1]

result = np.zeros((len(cell_positions), len(genes_interes)), dtype=np.float32)
gene_map = {g: j for j, g in enumerate(gene_positions)}

print("Leyendo datos...")
for i, (s, e) in enumerate(zip(all_starts, all_ends)):
    if e > s:
        col_idx = f['X']['indices'][s:e]
        values = f['X']['data'][s:e]
        mask = np.isin(col_idx, gene_positions)
        if mask.any():
            for c, v in zip(col_idx[mask], values[mask]):
                result[i, gene_map[c]] = v
    if i % 50000 == 0:
        print(f"  {i}/{len(cell_positions)}")

df_acc = pd.DataFrame(result, columns=genes_interes)
df_acc['subclass'] = meta_filtrada.loc[cell_names[cell_positions], 'subclass_short'].values
acc_mean = df_acc.groupby('subclass')[genes_interes].mean()
print(acc_mean)

# %%
df_dot = []
for subclass in df_acc['subclass'].unique():
    sub = df_acc[df_acc['subclass'] == subclass]
    for g in genes_interes:
        df_dot.append({
            'subclass': subclass,
            'gene': g,
            'mean': sub[g].mean(),
            'frac': (sub[g] > 0).mean()
        })

df_dot = pd.DataFrame(df_dot)

print(df_dot.groupby('gene')[['mean','frac']].describe())

df_rna = pd.read_csv("data/dot_plot_rnaseq_data.csv", index_col=0)

# %%
def make_dotplot_rna(df, gene_col, subclass_col, mean_col, frac_col,
                     title, output, midpoint=0.75, vmin=0, vmax=2.5,
                     frac_scale=50, figsize=(7,5)):
    
    gene_order = ['Rpl22', 'Rpl22l1', 'Rpl7', 'Rpl7l1', 'Rps27', 'Rps27l']
    subclass_order = sorted(df[subclass_col].unique())
    
    cmap_custom = LinearSegmentedColormap.from_list('blue_gray_red', ['blue', 'lightgray', 'red'])
    norm = TwoSlopeNorm(vmin=vmin, vcenter=midpoint, vmax=vmax)
    
    df = df.copy()
    df[frac_col] = (df[frac_col].clip(lower=70) - 70) / (100 - 70)
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_position([0.3, 0.15, 0.45, 0.75])
    
    for _, row in df.iterrows():
        x = gene_order.index(row[gene_col])
        y = subclass_order.index(row[subclass_col])
        ax.scatter(x=x, y=y,
                   s=row[frac_col] * frac_scale + 20,
                   c=[row[mean_col]],
                   cmap=cmap_custom, norm=norm,
                   edgecolors='gray', linewidths=0.3)
    
    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=45, ha='right', fontstyle='italic', fontsize=8)
    ax.set_yticks(range(len(subclass_order)))
    ax.set_yticklabels(subclass_order, fontsize=8)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    sm = plt.cm.ScalarMappable(cmap=cmap_custom, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.3, aspect=10, anchor=(0, 1.0))
    cbar.set_label(title, fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    
    for pct in [70, 80, 90, 100]:
        frac_norm = (pct - 70) / (100 - 70)
        ax.scatter([], [], s=frac_norm * frac_scale + 20, c='black', label=f'{pct}%')
    ax.legend(title='% expressing', bbox_to_anchor=(1.25, 0.0),
              loc='lower left', borderaxespad=0, fontsize=7, labelspacing=2,
              frameon=False)
    
    plt.savefig(output, format='svg')
    print(f"Guardado: {output}")


def make_dotplot_atac(df, gene_col, subclass_col, mean_col, frac_col,
                      title, output, frac_scale=50, max_frac=0.30, figsize=(7,5)):
    
    gene_order = ['Rpl22', 'Rpl22l1', 'Rpl7', 'Rpl7l1', 'Rps27', 'Rps27l']
    subclass_order = sorted(df[subclass_col].unique())
    
    cmap_custom = LinearSegmentedColormap.from_list('blue_gray_red', ['blue', 'lightgray', 'red'])
    midpoint = df[mean_col].median()
    norm = TwoSlopeNorm(vmin=df[mean_col].min(), vcenter=midpoint, vmax=df[mean_col].max())
    
    df = df.copy()
    df[frac_col] = (df[frac_col].clip(upper=max_frac)) / max_frac
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_position([0.3, 0.15, 0.45, 0.75])
    
    for _, row in df.iterrows():
        x = gene_order.index(row[gene_col])
        y = subclass_order.index(row[subclass_col])
        ax.scatter(x=x, y=y,
                   s=row[frac_col] * frac_scale + 20,
                   c=[row[mean_col]],
                   cmap=cmap_custom, norm=norm,
                   edgecolors='gray', linewidths=0.3)
    
    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=45, ha='right', fontstyle='italic', fontsize=8)
    ax.set_yticks(range(len(subclass_order)))
    ax.set_yticklabels(subclass_order, fontsize=8)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    sm = plt.cm.ScalarMappable(cmap=cmap_custom, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.3, aspect=10, anchor=(0, 1.0))
    cbar.set_label(title, fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    
    for pct in [1, 10, 20, 30]:
        frac_norm = (pct / 100) / max_frac
        ax.scatter([], [], s=frac_norm * frac_scale + 20, c='black', label=f'{pct}%')
    ax.legend(title='% accessible', bbox_to_anchor=(1.25, 0.0),
              loc='lower left', borderaxespad=0, fontsize=7, labelspacing=2,
              frameon=False)
    
    plt.savefig(output, format='svg')
    print(f"Guardado: {output}")


n_rna = df_rna['id'].nunique()
n_atac = df_dot['subclass'].nunique()

make_dotplot_rna(
    df=df_rna,
    gene_col='features.plot',
    subclass_col='id',
    mean_col='avg.exp',
    frac_col='pct.exp',
    title='Average Expression',
    output='Figuras/Figure4a.svg',
    figsize=(7, n_rna * 0.2)
)

make_dotplot_atac(
    df=df_dot,
    gene_col='gene',
    subclass_col='subclass',
    mean_col='mean',
    frac_col='frac',
    title='Mean accessibility (logCPM)',
    output='Figuras/Figure4b.svg',
    figsize=(7, n_atac * 0.2)
)


