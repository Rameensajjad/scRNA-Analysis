import numpy as np
import pandas as pd
import anndata as ad
import anndata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.sparse import csr_matrix

print("anndata version:", ad.__version__)

# ── PART 1 ───────────────────────────────────────────────
print("\n" + "="*60)
print("PART 1 – Building an AnnData from scratch")
print("="*60)

# 1.1
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)
print("\nInitial AnnData:\n", adata)
print("adata.X:\n", adata.X)

# 1.2
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
print("\nFirst 10 obs_names:", adata.obs_names[:10].tolist())

# 1.3
subset_view = adata[["Cell_1", "Cell_10"], ["Gene_5", "Gene_1900"]]
print("\nSubset view:", subset_view)

# 1.4
ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)
print("\nadata.obs head:\n", adata.obs.head())
print("Updated AnnData:", adata)

# 1.5
bdata = adata[adata.obs.cell_type == "B"]
print("\nB-cell subset:", bdata)

# 1.6
adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))
print("\nobsm keys:", list(adata.obsm.keys()))
print("AnnData after obsm/varm:", adata)

# 1.7
adata.uns["random"] = [1, 2, 3]
print("\nadata.uns:", dict(adata.uns))

# 1.8
adata.layers["log_transformed"] = np.log1p(adata.X)
print("\nAnnData with layers:", adata)

# 1.9
df_log = adata.to_df(layer="log_transformed")
print("\nlog_transformed DataFrame (first 3 rows, first 5 cols):\n", df_log.iloc[:3, :5])

# 1.10
adata.write("my_results.h5ad", compression="gzip")
print("\nSaved to my_results.h5ad")

# 1.11
obs_meta = pd.DataFrame(
    {
        "time_yr":         np.random.choice([0, 2, 4, 8], adata.n_obs),
        "subject_id":      np.random.choice(["subject 1", "subject 2", "subject 4", "subject 8"], adata.n_obs),
        "instrument_type": np.random.choice(["type a", "type b"], adata.n_obs),
        "site":            np.random.choice(["site x", "site y"], adata.n_obs),
    },
    index=adata.obs.index,
)
adata = ad.AnnData(adata.X, obs=obs_meta, var=adata.var)
print("\nAnnData with rich obs:", adata)

view_5x2 = adata[:5, ["Gene_1", "Gene_3"]]
print("View (5x2):", view_5x2)

adata_subset = adata[:5, ["Gene_1", "Gene_3"]].copy()
print("Copy from view:", adata_subset)

print("Before:", adata[:3, "Gene_1"].X.toarray().tolist())
adata[:3, "Gene_1"].X = [0, 0, 0]
print("After zeroing:", adata[:3, "Gene_1"].X.toarray().tolist())

adata_sub3 = adata[:3, ["Gene_1", "Gene_2"]]
adata_sub3.obs["foo"] = range(3)
print("After implicit copy:", adata_sub3)

print("Cells from time_yr 2 or 4:\n", adata[adata.obs.time_yr.isin([2, 4])].obs.head())

# 1.12
adata_backed = ad.io.read_h5ad("my_results.h5ad", backed="r")
print("\nIs backed:", adata_backed.isbacked)
print("Backed filename:", adata_backed.filename)
adata_backed.file.close()

# ── PART 2 ───────────────────────────────────────────────
print("\n" + "="*60)
print("PART 2 – Exploring a real PBMC dataset")
print("="*60)

# 2.1
adata = anndata.read_h5ad("pbmc3k_simulated.h5ad")
print("\nLoaded AnnData:\n", adata)

# 2.2
print("\nadata.X:\n", adata.X)
print("Non-zero values sample:", adata.X.data[:5])
print("Non-zero indices sample:", adata.X.indices[:5])
print("Sparsity:", round(adata.X.nnz / np.prod(adata.X.shape), 4))

# 2.3
print("\nLayers:", list(adata.layers.keys()))
print("raw dtype:", adata.layers["raw"].dtype)
print("raw sample:", adata.layers["raw"].data[:5])

adata.layers["counts_per_million"] = adata.layers["raw"].copy()
sc.pp.normalize_total(adata, target_sum=10**6, layer="counts_per_million")
print("Layers after CPM:", list(adata.layers.keys()))

genes_of_interest = ["CD8A", "CD4", "KLRB1"]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 2.5))
sc.pl.matrixplot(adata, groupby="louvain_cell_types", var_names=genes_of_interest,
                 layer="counts_per_million", ax=ax1, show=False)
ax1.set_title("CPM normalization")
sc.pl.matrixplot(adata, groupby="louvain_cell_types", var_names=genes_of_interest,
                 layer="raw", ax=ax2, show=False)
ax2.set_title("Raw counts")
plt.tight_layout()
plt.savefig("pbmc_layers_comparison.png", dpi=120, bbox_inches="tight")
plt.close()
print("Saved pbmc_layers_comparison.png")

# 2.4
print("\nobs columns:", adata.obs.keys().tolist())
print("obs head:\n", adata.obs.head())
n_b = sum(adata.obs["louvain_cell_types"] == "B cells")
print(f"Number of B cells: {n_b}")
adata.obs["is_low_quality"] = adata.obs["percent_mito"] > 0.03
print("obs after quality flag:\n", adata.obs.head())
print("var head:\n", adata.var.head())

# 2.5
print("\nobs_names (first 5):", adata.obs_names[:5].tolist())
print("var_names (first 5):", adata.var_names[:5].tolist())
adata.var_names = adata.var["gene_ids"]
print("var after Ensembl re-index:\n", adata.var.head(3))
adata.var_names = adata.var["gene_names"]
print("var restored:\n", adata.var.head(3))

# 2.6
adata_small = adata[:5, ["LYZ", "FOS", "MALAT1"]]
print("\nadata_small shape:", adata_small.shape)
print("X dense:\n", adata_small.X.toarray())
print("raw layer:\n", adata_small.layers["raw"].toarray())
print("obs:\n", adata_small.obs)
print("var:\n", adata_small.var)

adata_hq = adata[~adata.obs["is_low_quality"], :]
print("\nHigh-quality cells:", adata_hq.n_obs)
print("adata_hq obs head:\n", adata_hq.obs.head())

# 2.7
print("\nobsm keys:", list(adata.obsm.keys()))
for key in adata.obsm:
    print(f"  {key}: shape {adata.obsm[key].shape}")

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
b_flag = adata.obs["louvain_cell_types"] == "B cells"
scatter_kw = dict(s=3, linewidth=0, cmap="coolwarm")
for ax, key, title in zip(axes, ["X_pca", "X_tsne", "X_umap"], ["PCA", "t-SNE", "UMAP"]):
    ax.scatter(adata.obsm[key][:, 0], adata.obsm[key][:, 1],
               c=b_flag.astype(int), **scatter_kw)
    ax.set_title(title)
    ax.axis("off")
    ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("pbmc_embeddings.png", dpi=120, bbox_inches="tight")
plt.close()
print("Saved pbmc_embeddings.png")

# 2.8
print("\nobsp keys:", list(adata.obsp.keys()))
print("distances_all shape:", adata.obsp["distances_all"].shape)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].imshow(adata.obsp["distances_all"])
axes[0].set_title("Random order")
axes[0].set_xlabel("Euclidean distance in PCA space")
reorder = np.argsort(adata.obs["louvain_cell_types"])
axes[1].imshow(adata[reorder, :].obsp["distances_all"])
axes[1].set_title("Sorted by cell type")
axes[1].set_xlabel("Euclidean distance in PCA space")
plt.tight_layout()
plt.savefig("pbmc_distance_matrix.png", dpi=120, bbox_inches="tight")
plt.close()
print("Saved pbmc_distance_matrix.png")

# 2.9
print("\nuns keys:", list(adata.uns.keys()))
print("louvain params:", adata.uns["louvain"])
print("louvain colors:", adata.uns["louvain_colors"])
print("PCA variance (first 5):", adata.uns["pca"]["variance"][:5])

# 2.10
adata_view = adata[:5, 5:10]
print("\nView:", adata_view)
print("Before modification:\n", adata_view.X.toarray())
adata.X[0, 7] = 99
print("After modifying parent:\n", adata_view.X.toarray())

adata_view_copy = adata_view.copy()
print("Explicit copy:", adata_view_copy)

adata_view.obs["new_column"] = "Test"
print("After obs assignment (now real AnnData):", adata_view)

adata.X[0, 8] = 98
print("Parent X[0,5:10]:\n",  adata.X[:5, 5:10].toarray())
print("Former view X:\n",     adata_view.X.toarray())
print("Explicit copy X:\n",   adata_view_copy.X.toarray())

print("\nDone.")
