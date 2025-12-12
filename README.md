# RP_scRNAseq_Neurons_Code

This repository contains all code and scripts required to reproduce the analyses and generate the figures for the article:

"Single-cell RNA-seq Reveals Distinct Ribosomal Protein Expression Profiles Among Neuronal Subtypes".

Data and Code Availability

All single-cell RNA-seq datasets analyzed in this study are publicly available:

Yao et al. (2021) – 10x Genomics Chromium and Smart-seq2 datasets (NeMO: dat-jb2f34y)

Jin et al. (2025) – NeMO dataset (dat-61kfys3)

Hing et al. (2024) – GEO dataset (GSE240975)

The Smart-seq2 dataset was downloaded directly as a Seurat object (Seurat.ss.rda) from the Allen Institute Brain Map portal:
https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/mouse-whole-cortex-and-hippocampus-smart-seq

Required Input Files and Directory Structure

Before running any script, populate the data/ directory with the following files:

data/

- Yao et al. (2021) 10xGenomics Dataset count matrix. Chromium_expression_matrix.hdf5
- Yao et al. (2021) 10xGenomics Dataset Metadata. CTX_Hip_anno_10x.csv
- Yao et al. (2021) Smartseq2 Dataset Seurat Object. Seurat.ss.rda
- Yao et al. (2021) Smartseq2 Dataset Metadata. CTX_Hip_anno_SSv4.csv
- Ribosomal Protein Genes. RP.csv
- Jin et al. (2025) Aging Dataset count matrix. Mouse_Aging_10Xv3_counts_20241115.h5ad
- Jin et al. (2025) Aging Dataset Metadata. metadata_aging 
- Hing et al. (2024) Stress Vulnerability Dataset. Place inside data/Hing_2024_PFC_stress/original_data/ all raw matrices downloaded from GEO (GSE240975).
- Immunohistochemistry Quantification. Place inside data/Immuno_Quantification/ the three immunofluorescence quantification tables. 
    - Cuantificacion_Rps27_NeuN_Gad67.xlsx
    - Cuantificacion_Rpl29_NeuN_Gad67.xlsx
    - Cuantificacion_Rpl21_NeuN_Gad67.xlsx

The project requires a specific execution order.
The figure-generation scripts depend on intermediate files produced during preprocessing and differential expression analysis.
The complete pipeline is organized into:

- Dataset preprocessing
- Stability (SCSeg) analysis
- Differential expression analysis
- Figure generation
