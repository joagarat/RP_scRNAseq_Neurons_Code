#Smartseq

library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(RColorBrewer)
library(svglite)
library(ggrepel)
library(purrr)
library(patchwork)


load("data/Smartseq2dataset_seurat_filtered.RData")
RP = read.csv("data/RP.csv")
RP = as.character(RP$x)
RP =gsub("Rack1", "Gnb2l1", x = RP)
seurat_merged_RP = seurat_merged[RP, ]
seurat_merged_RP = SCTransform(seurat_merged_RP, vst.flavor = "v1")
seurat_merged_RP = RunPCA(seurat_merged_RP)
seurat_merged_RP = RunUMAP(seurat_merged_RP, dims = 1:8)
seurat_merged_RP = FindNeighbors(seurat_merged_RP, dims = 1:8)
seurat_merged_RP = FindClusters(seurat_merged_RP, resolution = 0.2)
ElbowPlot(seurat_merged_RP)
p_elbow <- ElbowPlot(seurat_merged_RP)
p_elbow_styled <- p_elbow +
  labs(
    title = NULL,                     
    x = "Principal Component",     
    y = "Standard Deviation"          
  ) +
  theme_classic() 

print(p_elbow_styled)

dir.create("Figuras/Figura_S2", recursive = TRUE, showWarnings = FALSE)
ggsave(
  plot = p_elbow_styled,
  filename = "Figuras/Figura_S2/FigS2_A.svg",
  width = 2.5,
  height = 2.5, 
  units = "in"
)

cell_counts_per_cluster <- table(seurat_merged_RP$SCT_snn_res.0.2)
clusters_with_at_least_10_cells <- names(cell_counts_per_cluster[cell_counts_per_cluster >= 10])
seurat_merged_RP_subset <- subset(seurat_merged_RP, idents = clusters_with_at_least_10_cells)
p_class <- DimPlot(seurat_merged_RP_subset, group.by = "class_label", raster = TRUE) +
  labs(
    title = NULL,         
    color = "Neuron Class"    
  ) +
  guides(color = guide_legend(
    title.position = "top", 
    ncol = 2,               
    override.aes = list(size = 2) 
  )) +
  theme_classic() +
  theme(
    legend.position = "bottom",   
    legend.title.align = 0.5,   
    
    aspect.ratio = 1, 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"), 
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(t = 0, b = 0) 
  )

print(p_class)
ggsave(
  filename = "Figuras/Figura_S2/FigS2_B.svg",
  plot = p_class,
  width = 1.9,  
  height = 2.75, 
  units = "in"
)

table_cluster_class <- table(seurat_merged_RP_subset$SCT_snn_res.0.2, seurat_merged_RP_subset$class_label)
df_cluster_class <- as.data.frame(table_cluster_class)
colnames(df_cluster_class) <- c("Cluster", "ClassLabel", "Count")
majority_labels_df <- df_cluster_class %>%
  group_by(Cluster) %>%
  mutate(TotalCount = sum(Count)) %>%
  ungroup() %>%
  mutate(Percentage = (Count / TotalCount) * 100) %>%
  group_by(Cluster) %>%
  slice_max(order_by = Count, n = 1) %>%
  ungroup() %>%
  mutate(LegendLabel = sprintf("Cluster %s (%d%% %s)", Cluster, round(Percentage), ClassLabel))
label_map <- setNames(majority_labels_df$LegendLabel, majority_labels_df$Cluster)
color_gaba <- "#F8766D" 
color_glut1 <- "#00BFC4"
color_glut2 <- "#008080" 

glut_clusters <- majority_labels_df %>%
  filter(ClassLabel == "Glutamatergic") %>%
  pull(Cluster)
color_map <- setNames(
  sapply(levels(seurat_merged_RP_subset$SCT_snn_res.0.2), function(cluster_level) {
    major_class <- majority_labels_df$ClassLabel[majority_labels_df$Cluster == cluster_level]
    if (major_class == "GABAergic") {
      return(color_gaba)
    } else if (cluster_level == glut_clusters[1]) {
      return(color_glut1)
    } else {
      return(color_glut2)
    }
  }),
  levels(seurat_merged_RP_subset$SCT_snn_res.0.2)
)
p_cluster <- DimPlot(seurat_merged_RP_subset, group.by = "SCT_snn_res.0.2", raster = TRUE) +
  scale_color_manual(values = color_map, labels = label_map) +
  labs(
    title = NULL,
    color = "Cluster" 
  ) +
  guides(color = guide_legend(
    title.position = "top", 
    ncol = 1,            
    override.aes = list(size = 3) 
  )) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title.align = 0.5,
    aspect.ratio = 1, # Fuerza 1:1
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"), 
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(t = 0, b = 0) 
  )

print(p_cluster)
dir.create("Figuras/Figura_2", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = "Figuras/Figura_2/Fig_2A.svg",
  plot = p_cluster,
  width = 2.5,
  height = 5,
  units = "in"
)

p_subclass <- DimPlot(seurat_merged_RP_subset, group.by = "subclass_label", raster = TRUE) +
  labs(
    title = NULL,
    color = "Neuron Subclass" 
  ) +
  guides(color = guide_legend(
    title.position = "top", 
    ncol = 3, 
    override.aes = list(size = 2) 
  )) +
  theme_classic() +
  theme(
    legend.position = "bottom", 
    legend.title.align = 0.5, 
    aspect.ratio = 1, 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"), 
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(t = 0, b = 0) 
  )

print(p_subclass)

ggsave(
  filename = "Figuras/Figura_S2/FigS2_C.svg",
  plot = p_subclass,
  width = 2.5,  
  height = 3.5, 
  units = "in"
)

regions <- c("ACA", "AI", "ALM", "HIP", "MOp", "RSP", "RSPv", "SSp", "SSs", "VIS", "VISp")
base_dir <- "data/SCSegIdx_smartseq_All_Regions/"
file_paths <- paste0(base_dir, regions, "/scSEG.csv")
segIdx_loaded <- map(file_paths, read.csv)
names(segIdx_loaded) <- regions
RP_df <- read.csv("data/RP.csv")
RP <- as.character(RP_df$x)
RP <- gsub("Rack1", "Gnb2l1", x = RP)
combined_data_raw <- bind_rows(segIdx_loaded, .id = "region")
combined_data <- combined_data_raw %>%
  mutate(gene_type = ifelse(toupper(gene) %in% toupper(RP), "RP", "Other")) %>%
  mutate(f_stats = pmax(0, f_stats))
gene_summary <- combined_data %>%
  group_by(gene, gene_type) %>%
  summarise(
    median_segIdx = median(segIdx, na.rm = TRUE),
    median_fstat = median(f_stats, na.rm = TRUE)
  ) %>%
  ungroup()
gene_summary <- gene_summary %>%
  mutate(gene_type = recode(gene_type, "Other" = "All Genes"))
color_palette <- c("RP" = "#d62728", "All Genes" = "gray50")
plot_A_density_segIdx <- ggplot(gene_summary, aes(x = median_segIdx, fill = gene_type)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(
    values = color_palette,
    name = "Gene Type" 
  ) +
  scale_x_reverse(limits = c(1, 0)) +
  labs(
    title = NULL, 
    x = "SegIdx", 
    y = "Density of Genes",
    caption = "Stable → Variable" 
  ) +
  theme_classic() +
  theme(
    legend.position = "top", # Mueve la leyenda arriba
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  )

print(plot_A_density_segIdx)


ggsave(
  filename = "Figuras/Figura_2/Fig_2B.svg",
  plot = plot_A_density_segIdx,
  width = 2.3,  # <-- Ancho de 1/3 de página
  height = 2.3, # <-- Alto IDÉNTICO al Panel A
  units = "in"
)

data_ranking_segIdx <- gene_summary %>%
  filter(!is.na(median_segIdx)) %>%
  arrange(desc(median_segIdx)) %>% 
  mutate(rank = row_number())

labels_b <- data_ranking_segIdx %>%
  filter(gene_type == "RP") %>%
  slice(c(1:5, (n()-4):n()))

gene_summary <- gene_summary %>%
  mutate(gene_type = recode(gene_type, "Other" = "All Genes"))


stable_genes_list <- gene_summary %>%
  filter(gene_type == "RP") %>%
  arrange(desc(median_segIdx)) %>%
  slice(1:5) %>%
  pull(gene) %>%
  paste(collapse = "\n") 

stable_annotation_text <- paste("Most Stable RPs:", stable_genes_list, sep = "\n")


variable_genes_list <- gene_summary %>%
  filter(gene_type == "RP") %>%
  arrange(median_segIdx) %>%
  slice(1:5) %>%
  pull(gene) %>%
  paste(collapse = "\n") 

variable_annotation_text <- paste("Most Variable RPs:", variable_genes_list, sep = "\n")
library(dplyr) # Necesario para filtrar y ordenar
top_5_estables <- data_ranking_segIdx %>%
  filter(gene_type == "RP") %>%
  arrange(desc(median_segIdx)) %>%
  slice(1:5)
top_5_variables <- data_ranking_segIdx %>%
  filter(gene_type == "RP") %>%
  arrange(median_segIdx) %>%
  slice(1:5)
top_10_labels <- rbind(top_5_estables, top_5_variables)
plot_B_ranking_segIdx_final <- ggplot() +
  geom_point(data = data_ranking_segIdx, aes(x = median_segIdx, y = rank), color = "gray80", size = 0.8) +
  geom_point(data = filter(data_ranking_segIdx, gene_type == "RP"), aes(x = median_segIdx, y = rank), color = color_palette["RP"], size = 2) +
  geom_point(
    data = top_10_labels,
    aes(x = median_segIdx, y = rank),
    shape = 21, size = 3, fill = color_palette["RP"], color = "black", stroke = 0.8
  ) +
  geom_point(data = data_ranking_segIdx, aes(x = median_segIdx, y = rank, fill = gene_type), alpha = 0, shape = 21) +
  scale_x_reverse(limits = c(1, 0)) +
  scale_fill_manual(name = "Gene Type", values = color_palette, labels = c("All Genes", "RP")) +
  guides(fill = guide_legend(override.aes = list(
    alpha = 1,
    shape = 22,    
    color = "black",
    size = 3.5,      
    stroke = 0.5   
  ))) +
  labs(
    title = NULL,
    x = "SegIdx",
    y = "Gene Rank",
    caption = "Stable → Variable"
  ) +
  theme_classic() +
  theme(
    legend.position = "top", 
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  )

print(plot_B_ranking_segIdx_final)
ggsave(
  filename = "Figuras/Figura_2/Fig_2C.svg",
  plot = plot_B_ranking_segIdx_final,
  width = 2.3,
  height = 2.3,
  units = "in"
)
print(top_5_estables)
print(top_5_variables)

color_palette <- c("RP" = "#d62728", "All Genes" = "gray50")

plot_C_density_FStat <- ggplot(gene_summary, aes(x = median_fstat, fill = gene_type)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(
    values = color_palette,
    name = "Gene Type"
  ) + 
  labs(
    title = NULL, 
    x = "F-Stat", 
    y = "Density of Genes",
    caption = "Non-Specific → Subclass Specific" 
  ) +
  theme_classic() + 
  theme(
    legend.position = "top", 
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  ) +
  coord_cartesian(xlim = c(0, 11)) 
print(plot_C_density_FStat)
ggsave(
  filename = "Figuras/Figura_2/Fig_2D.svg",
  plot = plot_C_density_FStat,
  width = 2.3,  # <-- Idéntico a los otros
  height = 2.3, # <-- Idéntico a los otros
  units = "in"
)

data_ranking_fstat <- gene_summary %>%
  filter(!is.na(median_fstat)) %>%
  arrange(median_fstat) %>% 
  mutate(rank = row_number())

labels_b <- data_ranking_fstat %>%
  filter(gene_type == "RP") %>%
  slice(c(1:5, (n()-4):n()))

gene_summary <- gene_summary %>%
  mutate(gene_type = recode(gene_type, "Other" = "All Genes"))

variable_genes_list <- gene_summary %>%
  filter(gene_type == "RP") %>%
  arrange(desc(median_fstat)) %>%
  slice(1:5) %>%
  pull(gene) %>%
  paste(collapse = "\n") 
variable_annotation_text <- paste("Most Specific RPs:", variable_genes_list, sep = "\n")
stable_genes_list <- gene_summary %>%
  filter(gene_type == "RP") %>%
  arrange(median_fstat) %>%
  slice(1:5) %>%
  pull(gene) %>%
  paste(collapse = "\n") 

stable_genes_list <- gsub("Gnb2l1", "Rack1", stable_genes_list)
stable_annotation_text <- paste("Non-Specific RPs:", stable_genes_list, sep = "\n")
library(dplyr)
top_5_estables_fstat <- data_ranking_fstat %>%
  filter(gene_type == "RP") %>%
  arrange(median_fstat) %>% 
  slice(1:5)
top_5_variables_fstat <- data_ranking_fstat %>%
  filter(gene_type == "RP") %>%
  arrange(desc(median_fstat)) %>%
  slice(1:5)
top_10_labels_fstat <- rbind(top_5_estables_fstat, top_5_variables_fstat)
plot_D_ranking_fstat_final <- ggplot() +
  geom_point(data = data_ranking_fstat, aes(x = median_fstat, y = rank), color = "gray80", size = 0.8) +
  geom_point(data = filter(data_ranking_fstat, gene_type == "RP"), aes(x = median_fstat, y = rank), color = color_palette["RP"], size = 2) +
  geom_point(
    data = top_10_labels_fstat,
    aes(x = median_fstat, y = rank),
    shape = 21, size = 3, fill = color_palette["RP"], color = "black", stroke = 0.8
  ) +
  geom_tile(data = data_ranking_fstat, aes(x = median_fstat, y = rank, fill = gene_type), alpha = 0) +
  scale_fill_manual(name = "Gene Type", values = color_palette, labels = c("All Genes", "RP")) +
  labs(
    title = NULL,
    x = "F-Stat",
    y = "Gene Rank",
    caption = "Non-Specific → Subclass Specific" 
  ) +
  coord_cartesian(xlim = c(0, 11)) +
  theme_classic() +
  theme(
    legend.position = "top", 
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"), 
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  )
print(plot_D_ranking_fstat_final)
ggsave(
  filename = "Figuras/Figura_2/Fig_2E.svg",
  plot = plot_D_ranking_fstat_final,
  width = 2.3,
  height = 2.3,
  units = "in"
)
print(top_5_variables_fstat)
print(top_5_estables_fstat)
fstats_median_rp <- gene_summary %>%
  filter(gene_type == "RP", !is.na(median_fstat)) %>%
  pull(median_fstat)
fstats_median_no_rp <- gene_summary %>%
  filter(gene_type == "All Genes", !is.na(median_fstat)) %>%
  pull(median_fstat)
test_global_fstat <- wilcox.test(fstats_median_rp, fstats_median_no_rp, alternative = "greater")
print(test_global_fstat)
lista_regiones <- unique(combined_data$region)
resultados_wilcoxon_por_region <- data.frame(
  region = character(),
  p_value = double(),
  statistic_W = double()
)
for (reg in lista_regiones) {
  data_region <- combined_data %>%
    filter(region == reg, !is.na(f_stats))
  fstats_rp    <- data_region %>% filter(gene_type == "RP") %>% pull(f_stats)
  fstats_no_rp <- data_region %>% filter(gene_type == "Other") %>% pull(f_stats)
  if (length(fstats_rp) > 0 && length(fstats_no_rp) > 0) {
    test_regional <- wilcox.test(fstats_rp, fstats_no_rp, alternative = "greater")
    resultados_wilcoxon_por_region <- rbind(resultados_wilcoxon_por_region, data.frame(
      region = reg,
      p_value = test_regional$p.value,
      statistic_W = test_regional$statistic
    ))
  }
}

print(resultados_wilcoxon_por_region %>% arrange(p_value))

regions <- c("ACA", "AI", "ALM", "HIP", "MOp", "RSP", "RSPv", "SSp", "SSs", "VIS", "VISp")
base_dir <- "data/SCSegIdx_smartseq_All_Regions/"
file_paths <- paste0(base_dir, regions, "/scSEG.csv")
segIdx_loaded <- map(file_paths, read.csv)
names(segIdx_loaded) <- regions
RP_df <- read.csv("data/RP.csv")
RP <- as.character(RP_df$x)
segIdx_loaded <- map(segIdx_loaded, ~ .x %>% 
                       mutate(gene = gsub("Gnb2l1", "Rack1", gene)))

generate_all_plots_for_region <- function(region_data, region_name, base_dir) {
  
  message(paste("--- Procesando región:", region_name, "---"))
  
  # --- A. Preparación de datos de la región ---
  gene_summary <- region_data %>%
    mutate(gene_type = ifelse(toupper(gene) %in% toupper(RP), "RP", "All Genes")) %>%
    mutate(f_stats = pmax(0, f_stats))
  
  color_palette <- c("RP" = "#d62728", "All Genes" = "gray50")
  
  # ==========================================================
  #         GRÁFICOS DE ESTABILIDAD (SegIdx)
  # ==========================================================
  
  # --- Plot A: Densidad SegIdx ---
  plot_A <- ggplot(gene_summary, aes(x = segIdx, fill = gene_type)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = color_palette, name = "Gene Type") +
    scale_x_reverse(limits = c(1, 0)) +
    labs(
      title = NULL, # <-- Título eliminado (se pone en Inkscape)
      x = "SegIdx", 
      y = "Density", # <-- Y-axis acortado
      caption = "Stable → Variable"
    ) +
    theme_classic() + # <-- base_size eliminado
    theme(
      legend.position = "none", # <-- Leyenda eliminada
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      plot.caption = element_text(hjust = 0.5, size = 8, face = "italic")
    )
  
  ggsave(
    filename = paste0("Figuras/Figura_S2/SegIdx_Density_", region_name, ".svg"), 
    plot = plot_A, 
    width = 2.5,  # <-- Ancho de 1/3 de página (7.5 / 3)
    height = 2.5  # <-- Alto = Ancho (para que sea cuadrado)
  )

  plot_C <- ggplot(gene_summary, aes(x = f_stats, fill = gene_type)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = color_palette, name = "Gene Type") +
    coord_cartesian(xlim = c(0, 11)) +
    labs(
      title = NULL, 
      x = "F-Stat", 
      y = "Density",
      caption = "Non-Specific → Subclass Specific" 
    ) +
    theme_classic() +
    theme(
      legend.position = "none", 
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      plot.caption = element_text(hjust = 0.5, size = 8, face = "italic")
    )
  ggsave(
    filename = paste0("Figuras/Figura_S2/FStat_Density_", region_name, ".svg"), 
    plot = plot_C, 
    width = 2.5,  
    height = 2.5  
  )
  
}

iwalk(segIdx_loaded, ~ generate_all_plots_for_region(.x, .y, base_dir = base_dir))

generate_all_plots_for_region <- function(region_data, region_name, base_dir) {
  color_palette <- c("RP" = "#d62728", "All Genes" = "gray50")
  data_ranking_segIdx <- region_data %>% 
    mutate(gene_type = ifelse(toupper(gene) %in% toupper(RP), "RP", "All Genes")) %>%
    arrange(desc(segIdx)) %>% 
    mutate(rank = row_number())
  stable_genes <- head(filter(data_ranking_segIdx, gene_type == "RP") %>% arrange(desc(segIdx)), 5)
  variable_genes <- head(filter(data_ranking_segIdx, gene_type == "RP") %>% arrange(segIdx), 5)
  stable_text <- paste("Most Stable:", paste(stable_genes$gene, collapse = "\n"), sep = "\n")
  variable_text <- paste("Most Variable:", paste(variable_genes$gene, collapse = "\n"), sep = "\n")
  full_segidx_text <- paste(stable_text, variable_text, sep = "\n\n")
  ranking_plot_segidx <- ggplot() +
    geom_point(data = data_ranking_segIdx, aes(x = segIdx, y = rank), color = "gray80", size = 0.8) +
    geom_point(data = filter(data_ranking_segIdx, gene_type == "RP"), aes(x = segIdx, y = rank), color = "#d62728", size = 2) +
    scale_x_reverse(limits = c(1, 0)) +
    coord_cartesian(ylim = c(0, 15000)) + 
    labs(
      title = NULL, x = "SegIdx", y = "Gene Rank", caption = "Stable → Variable"
    ) +
    theme_classic() + 
    theme(
      legend.position = "none", 
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"),
      plot.margin = margin(t = 5, r = 0, b = 5, l = 5, unit = "pt")
    )
  text_plot_segidx <- ggplot() +
    annotate("text", x = 0, y = 1, label = full_segidx_text, 
             hjust = 0, vjust = 1,
             size = 2.5) + 
    theme_void() +
    ylim(c(0, 1)) +
    theme(plot.margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt"))
  final_plot_B <- ranking_plot_segidx + text_plot_segidx + plot_layout(widths = c(1.5, 1)) 
  ggsave(
    filename = paste0("Figuras/Figura_S2/SegIdx_RankPlot_", region_name, ".svg"), 
    plot = final_plot_B, 
    width = 4,   
    height = 2.7  
  )
  data_ranking_fstat <- region_data %>% 
    mutate(gene_type = ifelse(toupper(gene) %in% toupper(RP), "RP", "All Genes")) %>%
    arrange(f_stats) %>% 
    mutate(rank = row_number())
  specific_genes <- head(filter(data_ranking_fstat, gene_type == "RP") %>% arrange(desc(f_stats)), 5)
  nonspecific_genes <- head(filter(data_ranking_fstat, gene_type == "RP") %>% arrange(f_stats), 5)
  specific_text <- paste("Most Specific:", paste(specific_genes$gene, collapse = "\n"), sep = "\n")
  nonspecific_text <- paste("Non-Specific:", paste(nonspecific_genes$gene, collapse = "\n"), sep = "\n")
  full_fstat_text <- paste(specific_text, nonspecific_text, sep = "\n\n")
  ranking_plot_fstat <- ggplot() +
    geom_point(data = data_ranking_fstat, aes(x = f_stats, y = rank), color = "gray80", size = 0.8) +
    geom_point(data = filter(data_ranking_fstat, gene_type == "RP"), aes(x = f_stats, y = rank), color = "#d62728", size = 2) +
    coord_cartesian(xlim = c(0, 11), ylim = c(0, 15000)) + 
    labs(
      title = NULL, x = "F-Stat", y = "Gene Rank", caption = "Non-Specific → Subclass Specific" 
    ) +
    theme_classic() + 
    theme(
      legend.position = "none", 
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"),
      plot.margin = margin(t = 5, r = 0, b = 5, l = 5, unit = "pt")
    )
  text_plot_fstat <- ggplot() +
    annotate("text", x = 0, y = 1, label = full_fstat_text, 
             hjust = 0, vjust = 1, 
             size = 2.5) + 
    theme_void() +
    ylim(c(0, 1)) +
    theme(plot.margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt"))
  final_plot_D <- ranking_plot_fstat + text_plot_fstat + plot_layout(widths = c(1.5, 1))
  
  ggsave(
    filename = paste0("Figuras/Figura_S2/FStat_RankPlot_", region_name, ".svg"), 
    plot = final_plot_D, 
    width = 4,    
    height = 2.7  
  )
}

iwalk(segIdx_loaded, ~ generate_all_plots_for_region(.x, .y, base_dir = base_dir))
