import scanpy as sc

adata = sc.read_h5ad("data/Mouse_Aging_10Xv3_counts_20241115.h5ad")

broad_roi_filter = adata.obs["broad_roi"].isin(["Isocortex"])
class_label_filter = ~adata.obs["class_label"].isin(["Astro-Epen", "Immune", "Vascular"])

adata_filtered = adata[broad_roi_filter & class_label_filter]

adata_filtered.write_h5ad("data/ArchivoFiltrado_Isocortex.h5ad")
