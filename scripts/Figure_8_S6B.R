library(ggplot2)
library(patchwork)
library(purrr)
library(data.table)
library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(svglite)
library(sigclust2)

RP <- read.csv("data/RP.csv")
RP <- as.character(RP$x)
RP <- gsub("Rack1", "Gnb2l1", x = RP) # Estandarizamos Rack1
path_stress <- "data/Hing_2024_PFC_stress_DEGs_ALL/"
filelist_stress <- list.files(path = path_stress, pattern="*.csv$")
df_list_stress <- lapply(filelist_stress, function(f) fread(file.path(path_stress, f)))
names(df_list_stress) <- gsub(filelist_stress, pattern="\\..*", replacement="")
df_stress <- bind_rows(df_list_stress, .id = "id")
df_stress_clean <- na.omit(df_stress)

path_control <- "data/Chromium_DEGs_PL_ILA_ORB/PL-ILA-ORB/"
filelist_control <- list.files(path = path_control, pattern="*.csv$")
df_list_control <- lapply(filelist_control, function(f) fread(file.path(path_control, f)))
names(df_list_control) <- gsub(filelist_control, pattern="\\..*", replacement="")
df_control <- bind_rows(df_list_control, .id = "id")
df_control_clean <- na.omit(df_control)
df_control_clean$id <- gsub("_vs_", "vs", df_control_clean$id)

ordenar_id <- function(id) {
  subtipos <- unlist(strsplit(id, "vs"))
  id_ordenado <- paste(sort(subtipos), collapse = "vs")
  invertido <- ifelse(id == id_ordenado, FALSE, TRUE)
  return(c(id_ordenado, invertido))
}

df_control_std <- df_control_clean %>%
  mutate(id_info = map(id, ordenar_id)) %>%
  mutate(id_estandarizado = sapply(id_info, `[`, 1),
         invertido = as.logical(sapply(id_info, `[`, 2))) %>%
  select(-id_info)

df_stress_std <- df_stress_clean %>%
  mutate(id_info = map(id, ordenar_id)) %>%
  mutate(id_estandarizado = sapply(id_info, `[`, 1),
         invertido = as.logical(sapply(id_info, `[`, 2))) %>%
  select(-id_info)

tabla_comun <- inner_join(
  df_control_std, 
  df_stress_std, 
  by = c("id_estandarizado", "V1"),
  suffix = c("_control", "_stress")
)



tabla_comun <- tabla_comun %>%
  mutate(
    log2FoldChange_control = ifelse(invertido_control, -log2FoldChange_control, log2FoldChange_control),
    log2FoldChange_stress = ifelse(invertido_stress, -log2FoldChange_stress, log2FoldChange_stress)
  )

tabla_comun <- tabla_comun %>%
  mutate(gene_type = ifelse(toupper(V1) %in% toupper(RP), "RP", "Non-RP"))

max_abs_logfc <- max(abs(c(tabla_comun$log2FoldChange_control, tabla_comun$log2FoldChange_stress)), na.rm = TRUE) * 1.05

correlacion_global <- cor(tabla_comun$log2FoldChange_control, tabla_comun$log2FoldChange_stress, method = "spearman")

tabla_comun_plot <- tabla_comun %>%
  mutate(gene_type = ifelse(toupper(V1) %in% toupper(RP), "RP", "Non-RP")) %>%
  arrange(gene_type)

max_abs_logfc <- max(abs(c(tabla_comun_plot$log2FoldChange_control, tabla_comun_plot$log2FoldChange_stress)), na.rm = TRUE) * 1.05

tabla_comun_plot <- tabla_comun %>%
  mutate(gene_type = ifelse(toupper(V1) %in% toupper(RP), "RP", "Non-RP")) %>%
  arrange(gene_type) 

logfc_threshold <- log2(1.5)
padj_threshold <- 0.05

scatter_data_rps_stress <- tabla_comun %>%
  filter(toupper(V1) %in% toupper(RP)) %>%
  mutate(
    color_category = case_when(
      padj_control < padj_threshold & abs(log2FoldChange_control) > logfc_threshold &
        padj_stress < padj_threshold & abs(log2FoldChange_stress) > logfc_threshold &
        sign(log2FoldChange_control) == sign(log2FoldChange_stress) ~ "Concordant RP",
      padj_control < padj_threshold & abs(log2FoldChange_control) > logfc_threshold &
        padj_stress < padj_threshold & abs(log2FoldChange_stress) > logfc_threshold &
        sign(log2FoldChange_control) != sign(log2FoldChange_stress) ~ "Discordant RP",
      TRUE ~ "Non-Replicated RP"
    ),
    color_category = factor(color_category, levels = c("Non-Replicated RP", "Discordant RP", "Concordant RP"))
  ) %>%
  arrange(color_category)

max_abs_logfc_rp <- max(abs(c(scatter_data_rps_stress$log2FoldChange_control, scatter_data_rps_stress$log2FoldChange_stress)), na.rm = TRUE) * 1.05
correlacion_rp <- cor(scatter_data_rps_stress$log2FoldChange_control, scatter_data_rps_stress$log2FoldChange_stress, method = "spearman")

library(patchwork)

panel_a_global_corr <- ggplot(tabla_comun_plot, aes(x = log2FoldChange_control, y = log2FoldChange_stress)) +
  geom_point(aes(color = gene_type), alpha = 0.6, size = 0.6) +
  scale_color_manual(name = "Gene Set", values = c("RP" = "#e41a1c", "Non-RP" = "grey80")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-max_abs_logfc, max_abs_logfc), ylim = c(-max_abs_logfc, max_abs_logfc)) +
  
  annotate("text", x = -max_abs_logfc * 0.95, y = max_abs_logfc * 0.95, 
           label = paste("Spearman's ρ =", round(correlacion_global, 2)), 
           hjust = 0, vjust = 1, size = 2.5) + 
  
  theme_classic() + 
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,-5,-5,-5),
    aspect.ratio = 1,
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(title = NULL, x = "Log2FC (Yao et al., 2021)", y = "Log2FC (Hing et al., 2024)")

panel_b_rp_corr_stress <- ggplot(scatter_data_rps_stress, aes(x = log2FoldChange_control, y = log2FoldChange_stress)) +
  geom_point(aes(color = color_category), size = 0.6, alpha = 0.8) +
  scale_color_manual(name = "RP Concordance", 
                     values = c("Concordant RP" = "#e41a1c", "Discordant RP" = "#377eb8", "Non-Replicated RP" = "grey80")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-max_abs_logfc_rp, max_abs_logfc_rp), ylim = c(-max_abs_logfc_rp, max_abs_logfc_rp)) +
  
  annotate("text", x = -max_abs_logfc_rp * 0.95, y = max_abs_logfc_rp * 0.95, 
           label = paste("Spearman's ρ =", round(correlacion_rp, 2)), 
           hjust = 0, vjust = 1, size = 2.5) + # Tamaño estandarizado IDÉNTICO al Panel A
  
  theme_classic() + 
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,-5,-5,-5),
    aspect.ratio = 1,
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7)
  ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(alpha = 1, size = 3))) +
  labs(title = NULL, x = "Log2FC (Yao et al., 2021)", y = NULL)
dir.create("Figuras/Figura_8", recursive = TRUE, showWarnings = FALSE)
ggsave("Figuras/Figura_8/Figure_8A.png", plot = panel_a_global_corr, width = 2.5, height = 2.8, units = "in", dpi = 300)
ggsave("Figuras/Figura_8/Figure_8B.png", plot = panel_b_rp_corr_stress, width = 2.5, height = 2.8, units = "in", dpi = 300)

final_figure_AB <- panel_a_global_corr + panel_b_rp_corr_stress
ggsave("Figuras/Figura_8/Figure_8AB.png", plot = final_figure_AB, width = 4.0, height = 2.5, units = "in", dpi = 300)

library(ggrepel)

path_stress_effect <- "data/Hing_2024_PFC_stress_DEGs_High_vs_Low/"
filelist_stress_effect <- list.files(path = path_stress_effect, pattern="*.csv$")
df_list_stress_effect <- lapply(filelist_stress_effect, function(f) fread(file.path(path_stress_effect, f)))
df_stress_effect <- bind_rows(df_list_stress_effect, .id = "subclass_comparison")
df_stress_effect_clean <- na.omit(df_stress_effect)

volcano_data_stress <- df_stress_effect_clean %>%
  mutate(
    is_RP = toupper(V1) %in% toupper(RP),
    color_category = case_when(
      padj < 0.05 & abs(log2FoldChange) > 0.58 & is_RP ~ "RP DEGs",
      padj < 0.05 & abs(log2FoldChange) > 0.58 & !is_RP ~ "Non-RP DEGs",
      TRUE ~ "Not Significant"
    ),
    color_category = factor(color_category, levels = c("Not Significant", "Non-RP DEGs", "RP DEGs"))
  ) %>%
  arrange(color_category) 

genes_a_etiquetar <- volcano_data_stress %>% 
  filter(color_category == "RP DEGs") %>%
  arrange(padj) %>%
  slice_head(n = 10)

panel_c_volcano_stress <- ggplot(volcano_data_stress, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color_category), size = 0.8, alpha = 0.7) +
  scale_color_manual(
    name = NULL, 
    values = c(
      "RP DEGs" = "#e41a1c", 
      "Non-RP DEGs" = "#FCCDE5", 
      "Not Significant" = "grey80"
    )
  ) +
  geom_text_repel(
    data = genes_a_etiquetar, 
    aes(label = V1), 
    size = 2, 
    box.padding = 0.3, 
    point.padding = 0.1, 
    min.segment.length = 0, 
    max.overlaps = Inf,
    segment.size = 0.2
  ) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.3) +
  theme_classic() + 
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,-5,-5,-5),
    axis.title = element_text(size = 9), 
    axis.text = element_text(size = 7)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = NULL,
    x = "Log2FC (High/Low Stress Vulnerability)", 
    y = "-log10(p-adj)"
  )

print(panel_c_volcano_stress)

ggsave(
  filename = "Figuras/Figura_8/Figure_8C.png", 
  plot = panel_c_volcano_stress, 
  width = 3,  
  height = 2.5, 
  units = "in", 
  dpi = 300
)

seurat_stress <- LoadSeuratRds("data/Hing_2024_PFC_stress/SeuratProject.Rds")
RP <- read.csv("data/RP.csv")
RP <- as.character(RP$x)
inhibitory_subclasses <- c("Sst", "Pvalb", "Vip", "Lamp5", "Sncg", "Sst Chodl")
excitatory_subclasses <- c("L5 PT CTX", "L6 CT CTX", "L4/5 IT CTX", "CR","L5 NP CTX","L2/3 IT CTX","L6 IT CTX", "L5 IT TPE-ENT","L5 IT CTX","L6b CTX") 

all_subclasses_list <- unique(seurat_stress$predicted.id)
subclass_info <- data.frame(subclass = all_subclasses_list) %>%
  mutate(
    Class = ifelse(subclass %in% inhibitory_subclasses, "GABAergic", "Glutamatergic")
  )
rownames(subclass_info) <- subclass_info$subclass

seurat_stress$stress_group <- str_extract(seurat_stress$replica, "^[^_]+")
seurat_low <- subset(seurat_stress, subset = stress_group == "low")
seurat_medium <- subset(seurat_stress, subset = stress_group == "medium")
seurat_high <- subset(seurat_stress, subset = stress_group == "high")
calculate_avg_expression_matrix <- function(seurat_object, rp_genes) {
  avg_expr <- AverageExpression(
    object = seurat_object,
    features = rp_genes,
    group.by = "predicted.id",
    assay = "RNA"
  )
  return(avg_expr$RNA)
}
m_low <- calculate_avg_expression_matrix(seurat_low, RP)
m_medium <- calculate_avg_expression_matrix(seurat_medium, RP)
m_high <- calculate_avg_expression_matrix(seurat_high, RP)

common_genes <- intersect(rownames(m_low), intersect(rownames(m_medium), rownames(m_high)))
common_subtypes <- intersect(colnames(m_low), intersect(colnames(m_medium), colnames(m_high)))

m_low <- m_low[common_genes, common_subtypes]
m_medium <- m_medium[common_genes, common_subtypes]
m_high <- m_high[common_genes, common_subtypes]
colnames(m_low) <- paste0(colnames(m_low), "_low")
colnames(m_medium) <- paste0(colnames(m_medium), "_medium")
colnames(m_high) <- paste0(colnames(m_high), "_high")
m_combined <- cbind(m_low, m_medium, m_high)
m_combined_zscore_by_col <- scale(m_combined)
m_combined_zscore_by_col[is.na(m_combined_zscore_by_col)] <- 0

hc_cols <- hclust(dist(t(m_combined_zscore_by_col)), method = "ward.D2")
dend <- as.dendrogram(hc_cols)

k <- 3 
clusters <- cutree(hc_cols, k = k)
desired_cluster_order <- c(1, 2, 3) 
ordered_labels <- labels(dend)
clusters_in_dend_order <- clusters[ordered_labels] 
weights_map <- setNames(1:length(desired_cluster_order), desired_cluster_order)
reordering_weights <- weights_map[as.character(clusters_in_dend_order)]
dend_reordered <- reorder(dend, wts = reordering_weights, agglo.fun = "mean")

column_names <- colnames(m_combined_zscore_by_col)
annotation_df <- data.frame(sample = column_names, row.names = column_names) %>%
  separate(sample, into = c("Subclass", "Stress_temp"), sep = "_", remove = FALSE) %>% 
  mutate(
    Class = case_when(
      Subclass %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip") ~ "GABAergic",
      TRUE ~ "Glutamatergic"
    ),
    SigClust = as.numeric(factor(Subclass)) %% 8 + 1,
    `Stress Vulnerability` = Stress_temp 
  )

if (!exists("base_sigclust_colors")) {
  base_sigclust_colors <- RColorBrewer::brewer.pal(8, "Set2")
  names(base_sigclust_colors) <- 1:8
}

subclass_to_cluster_map <- annotation_df %>%
  dplyr::count(Subclass, SigClust) %>%
  group_by(Subclass) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup()

final_subclass_colors <- list()
for (cluster_id in unique(subclass_to_cluster_map$SigClust)) {
  base_color <- tryCatch(base_sigclust_colors[as.character(cluster_id)], error = function(e) "#CCCCCC")
  if(is.na(base_color)) base_color <- "#CCCCCC"
  
  subclasses_in_cluster <- subclass_to_cluster_map %>%
    filter(SigClust == cluster_id) %>%
    pull(Subclass)
  n_shades <- length(subclasses_in_cluster)
  color_generator <- colorRampPalette(c(adjustcolor(base_color, red.f = 1.2, green.f = 1.2, blue.f = 1.2), adjustcolor(base_color, red.f = 0.7, green.f = 0.7, blue.f = 0.7)))
  shades <- color_generator(n_shades)
  names(shades) <- subclasses_in_cluster
  final_subclass_colors[[length(final_subclass_colors) + 1]] <- shades
}
subclass_colors <- unlist(final_subclass_colors)
class_colors <- c("GABAergic" = "#882255", "Glutamatergic" = "#009988", "Unknown" = "grey")
stress_colors <- c("low" = "lightblue", "medium" = "orange", "high" = "firebrick")

annotation_colors <- list(
  Class = class_colors,
  Subclass = subclass_colors,
  `Stress Vulnerability` = stress_colors 
)

top_annotation <- HeatmapAnnotation(
  df = annotation_df %>% select(Class, Subclass), 
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  show_legend = TRUE,
  annotation_legend_param = list(
    Subclass = list(ncol = 3)
  )
)

bottom_annotation <- HeatmapAnnotation(
  df = annotation_df %>% select(`Stress Vulnerability`), 
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  show_legend = TRUE
)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

svglite("Figuras/Figura_8/Figura_8D_Heatmap.svg", 
        width = 10.5,  
        height = 9.5) 

ht <- Heatmap(m_combined_zscore_by_col,
              name = "Z-score",
              col = col_fun,
              top_annotation = top_annotation,
              bottom_annotation = bottom_annotation,
              cluster_columns = rev(dend_reordered), 
              cluster_rows = TRUE,
              show_column_names = FALSE, 
              row_names_gp = gpar(fontsize = 8),
              show_row_dend = TRUE,
              show_column_dend = TRUE
)

draw(ht, 
     heatmap_legend_side = "left",     
     annotation_legend_side = "bottom",
     merge_legend = TRUE 
)

dev.off()

final_column_order <- labels(dend_reordered)
m_reordered_for_shc <- m_combined_zscore_by_col[, final_column_order]
shc_result_reordered <- shc(as.matrix(t(m_reordered_for_shc)), metric = "euclidean", linkage = "ward.D2")
dir.create("Figuras/Figure_S6", recursive = TRUE, showWarnings = FALSE)
svglite(
  "Figuras/Figure_S6/Figure_S6B.svg", 
  width = 15, 
  height = 10
)

plot(shc_result_reordered, hang = -1, main = "Hierarchical Clustering with Significance (Reordered)")

dev.off()

