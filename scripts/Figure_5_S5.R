
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(tibble)
library(svglite)
library(ggbeeswarm)
library(RColorBrewer) 

DEG_dir = "data/Smartseq_DEGs/"
RP_df <- read.csv("data/RP.csv")
RP <- as.character(RP_df$x)
RP <- gsub("Rack1", "Gnb2l1", x = RP)
filelist <- list.files(path = DEG_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
regions <- sapply(stringr::str_split(filelist, "/"), function(path) path[length(path) - 1])
df_input_list <- lapply(seq_along(filelist), function(i) {
  df <- data.table::fread(filelist[i])
  df$region <- regions[i]
  df$comparison_raw <- gsub("^.*/|\\.csv$", "", filelist[i])
  df
})
df_allregions <- dplyr::bind_rows(df_input_list)
df_allregions_RP <- df_allregions[df_allregions$V1 %in% RP, ]
df_allregions_RP$log2FoldChange = df_allregions_RP$log2FoldChange * -1
standardize_comparison <- function(comparison_name, log2fc) {
  parts <- strsplit(comparison_name, "vs")[[1]]
  if (length(parts) == 2) {
    sorted_parts <- sort(parts)
    order_changed <- !identical(parts, sorted_parts)
    standardized_name <- paste(sorted_parts, collapse = "vs")
    standardized_log2fc <- if (order_changed) -log2fc else log2fc
    return(list(name = standardized_name, log2fc = standardized_log2fc))
  } else {
    return(list(name = comparison_name, log2fc = log2fc))
  }
}
df_corrected <- df_allregions_RP %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    standardized = list(standardize_comparison(comparison_raw, log2FoldChange)),
    comparison_id = standardized$name,
    log2FoldChange_std = standardized$log2fc
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(V1, region, comparison_id, log2FoldChange_std, padj)
inhibitory_subclasses <- c("Sst", "Pvalb", "Vip", "Lamp5", "Sncg", "Sst Chodl")
excitatory_subclasses <- c("CA1-ProS","CA3", "Car3", "DG","L2_3 IT CTX", "L2_3 IT PPP","L4 RSP-ACA",
                           "L4_5 IT CTX", "L5 IT CTX", "L5 NP CTX","L5 PT CTX","L6 CT CTX", 
                           "L6 IT CTX" , "L6b CTX")
all_subclasses_list <- c(inhibitory_subclasses, excitatory_subclasses)
subclass_info <- data.frame(subclass = all_subclasses_list) %>%
  mutate(Class = ifelse(subclass %in% inhibitory_subclasses, "Inhibitory", "Excitatory"))
rownames(subclass_info) <- subclass_info$subclass

top_genes <- df_corrected %>%
  filter(padj < 0.05, abs(log2FoldChange_std) > 0.58496250072) %>%
  count(V1, sort = TRUE) %>%
  slice_head(n = 25) %>%
  pull(V1)

winner_loser_data <- df_corrected %>%
  filter(padj < 0.05, abs(log2FoldChange_std) > 0.58496250072) %>%
  tidyr::separate(comparison_id, into = c("subclass1", "subclass2"), sep = "vs", remove = FALSE) %>%
  mutate(
    winner_subclass = case_when(
      log2FoldChange_std < 0 ~ subclass2,
      log2FoldChange_std > 0 ~ subclass1
    ),
    loser_subclass = case_when(
      log2FoldChange_std < 0 ~ subclass1,
      log2FoldChange_std > 0 ~ subclass2
    )
  )

total_comparisons_df <- data.frame(comparison_id = unique(df_corrected$comparison_id)) %>%
  tidyr::separate(comparison_id, into = c("subclass1", "subclass2"), sep = "vs", remove = FALSE) %>%
  tidyr::pivot_longer(cols = c("subclass1", "subclass2"), values_to = "subclass") %>%
  count(subclass, name = "total_comparisons")

win_counts <- winner_loser_data %>% group_by(V1, winner_subclass) %>% summarise(win_count = n(), .groups = 'drop') %>% rename(gene = V1, subclass = winner_subclass)
loss_counts <- winner_loser_data %>% group_by(V1, loser_subclass) %>% summarise(loss_count = n(), .groups = 'drop') %>% rename(gene = V1, subclass = loser_subclass)

win_normalized <- win_counts %>%
  left_join(total_comparisons_df, by = "subclass") %>%
  mutate(win_score = ifelse(total_comparisons > 0, win_count / total_comparisons, 0))

loss_normalized <- loss_counts %>%
  left_join(total_comparisons_df, by = "subclass") %>%
  mutate(loss_score = ifelse(total_comparisons > 0, loss_count / total_comparisons, 0))

all_genes <- top_genes
all_subclasses <- total_comparisons_df$subclass
win_complete <- win_normalized %>%
  filter(!is.na(subclass)) %>%
  tidyr::complete(gene = all_genes, subclass = all_subclasses, fill = list(win_score = 0))
loss_complete <- loss_normalized %>%
  filter(!is.na(subclass)) %>%
  tidyr::complete(gene = all_genes, subclass = all_subclasses, fill = list(loss_score = 0))

heatmap_b1_matrix_full <- win_complete %>%
  pivot_wider(id_cols = gene, names_from = subclass, values_from = win_score) %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.matrix()
heatmap_b1_matrix_full <- heatmap_b1_matrix_full[top_genes, ]

heatmap_b2_matrix_full <- loss_complete %>%
  pivot_wider(id_cols = gene, names_from = subclass, values_from = loss_score) %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.matrix()
heatmap_b2_matrix_full <- heatmap_b2_matrix_full[top_genes, ]

total_activity_matrix <- heatmap_b1_matrix_full + heatmap_b2_matrix_full
columnas_activas <- colSums(total_activity_matrix) > 0

subclases_removidas <- names(columnas_activas[!columnas_activas])
heatmap_b1_matrix <- heatmap_b1_matrix_full[, columnas_activas]
heatmap_b2_matrix <- heatmap_b2_matrix_full[, columnas_activas]

columnas_finales <- colnames(heatmap_b1_matrix)

annotation_df_prep <- subclass_info %>%
  filter(subclass %in% columnas_finales) %>%
  mutate(Class = recode(Class, 
                        "Inhibitory" = "GABAergic", 
                        "Excitatory" = "Glutamatergic")) %>%
  as.data.frame()
rownames(annotation_df_prep) <- annotation_df_prep$subclass

annotation_df_cols <- annotation_df_prep[, "Class", drop = FALSE]

colores_clase <- list(
  Class = c("GABAergic" = "#882255", "Glutamatergic" = "#009988")
)

annotation_df_cols <- annotation_df_cols[colnames(heatmap_b1_matrix), , drop = FALSE]

column_annotation <- HeatmapAnnotation(
  df = annotation_df_cols,
  col = colores_clase,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(0.5, "cm")
)

col_fun_wins <- colorRamp2(c(0, 0.1, 3), c("grey95", "#FEE0D2", "#99000D"))
col_fun_loss <- colorRamp2(c(0, 0.1, 3), c("grey95", "#DEEBF7", "#08306B"))

panel_b1_winners <- Heatmap(heatmap_b1_matrix,
                            name = "OE Score",
                            col = col_fun_wins,
                            top_annotation = column_annotation, 
                            cluster_rows = FALSE, 
                            cluster_columns = TRUE,
                            column_title = "Overexpression",
                            # --- AJUSTE DE FUENTE ---
                            row_names_gp = gpar(fontsize = 8), 
                            column_names_gp = gpar(fontsize = 8) 
)

pdf(NULL)
drawn_ht1 <- draw(panel_b1_winners)
master_column_order_names <- colnames(heatmap_b1_matrix)[column_order(drawn_ht1)]
dev.off()

reordering_weights_exp <- (1:length(master_column_order_names))^2
names(reordering_weights_exp) <- master_column_order_names
weights_for_heatmap2 <- reordering_weights_exp[colnames(heatmap_b2_matrix)]

panel_b2_losers <- Heatmap(heatmap_b2_matrix,
                           name = "UE Score",
                           col = col_fun_loss,
                           top_annotation = column_annotation, 
                           cluster_rows = FALSE, 
                           cluster_columns = TRUE,
                           column_dend_reorder = weights_for_heatmap2,
                           column_title = "Underexpression",
                           row_names_gp = gpar(fontsize = 8), 
                           column_names_gp = gpar(fontsize = 8)
)

panel_b_combined <- panel_b1_winners + panel_b2_losers
dir.create("Figuras/Figura_5", recursive = TRUE, showWarnings = FALSE)
svglite(
  file = "Figuras/Figura_5/Figura_5B.svg",
  width = 6.5,  
  height = 5.0 
)
draw(panel_b_combined)
dev.off()


top_3_genes <- df_corrected %>%
  filter(padj < 0.05, abs(log2FoldChange_std) > 0.58496250072) %>%
  count(V1, sort = TRUE) %>%
  slice_head(n = 3) %>%
  pull(V1)

plot_data_c <- df_corrected %>%
  filter(V1 %in% top_3_genes) %>%
  separate(comparison_id, into = c("subclass1", "subclass2"), sep = "vs", remove = FALSE) %>%
  mutate(
    type1 = subclass_info[subclass1, "Class"],
    type2 = subclass_info[subclass2, "Class"],
    `Comparison Type` = case_when(
      type1 == "Inhibitory" & type2 == "Inhibitory" ~ "Intra-Inhibitory",
      type1 == "Excitatory" & type2 == "Excitatory" ~ "Intra-Excitatory",
      TRUE ~ "Excitatory vs Inhibitory"
    ),
    abs_log2FC = abs(log2FoldChange_std),
    neg_log10_padj = -log10(padj),
    is_significant = ifelse(padj < 0.05 & abs_log2FC > 0.58496250072, "Significant", "Not Significant")
  ) %>%
  mutate(`Comparison Type` = factor(`Comparison Type`, 
                                    levels = c("Excitatory vs Inhibitory", "Intra-Excitatory", "Intra-Inhibitory")))
color_upper_limit <- plot_data_c %>%
  filter(is_significant == "Significant") %>%
  summarise(limit = quantile(neg_log10_padj, probs = 0.90, na.rm = TRUE)) %>%
  pull(limit)

if(is.na(color_upper_limit) || color_upper_limit < 10) color_upper_limit <- 50

panel_c_beeswarm_final <- ggplot(plot_data_c, aes(x = abs_log2FC, y = `Comparison Type`)) +
  
  geom_quasirandom(data = . %>% filter(is_significant == "Not Significant"),
                   color = "grey80", alpha = 0.5, size = 2) +
  
  geom_quasirandom(data = . %>% filter(is_significant == "Significant"),
                   aes(color = neg_log10_padj), size = 3) +
  
  facet_wrap(~ V1, nrow = 1) +
  
  geom_vline(xintercept = 0.58496250072, linetype = "dashed", color = "black") +
  
  scale_color_viridis_c(
    option = "plasma", 
    name = "-log10(p-adj)",
    limits = c(-log10(0.05), 20),
    oob = scales::squish 
  ) +
  
  theme_minimal() + 

  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(colour = "grey80", fill=NA),
    legend.position = "bottom",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 10) # <-- 10pt (no 14)
  ) +
  labs(
    title = NULL,
    x = "|log2 Fold Change|",
    y = NULL
  )
print(panel_c_beeswarm_final)

ggsave(
  "Figuras/Figura_5/Figure_5C.svg", 
  plot = panel_c_beeswarm_final, 
  width = 5.0,  
  height = 3.0 
)

df_allregions_RP_filtered = df_corrected[df_corrected$padj<0.05 & abs(df_corrected$log2FoldChange_std) > 0.58496250072, ]
region_counts <- df_allregions_RP_filtered %>%
  group_by(region) %>%
  summarise(
    total_rows = n(),
    unique_V1_count = n_distinct(V1)
  )
library(tidyr)
stats_RP <- df_allregions_RP_filtered %>%
  group_by(V1) %>%
  summarise(total_count = n())

region_counts_RP <- df_allregions_RP_filtered %>%
  group_by(V1, region) %>%
  summarise(count = n()) %>%
  spread(key = region, value = count, fill = 0)
final_stats <- left_join(stats_RP, region_counts_RP, by = "V1")


data_long <- final_stats %>%
  select(-total_count) %>%
  pivot_longer(
    cols = -V1,
    names_to = "region",
    values_to = "detections"
  )
colnames(data_long)[1] <- "gene"

gene_order <- final_stats %>%
  arrange(total_count) %>%
  pull(V1)

data_long_ordered <- data_long %>%
  mutate(gene = factor(gene, levels = gene_order)) 

numero_de_regiones <- n_distinct(data_long$region)
colores_nuevos <- brewer.pal(n = numero_de_regiones, name = "Set3")

panel_a_barplot <- ggplot(data_long_ordered, aes(x = detections, y = gene, fill = region)) +
  geom_col() +
  scale_fill_manual(values = colores_nuevos) +
  labs(
    title = NULL,
    x = "Total DE Findings",
    y = NULL,
    fill = "Brain Region"
  ) +
  theme_classic() + 
  theme(
    legend.position = "bottom",
    legend.title.position = "top", # <-- Pone el título de la leyenda arriba
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.text.y = element_text(face = "italic", size = 8) 
  )
print(panel_a_barplot)


ggsave(
  "Figuras/Figura_5/Figura_5A.svg", 
  plot = panel_a_barplot, 
  width = 4,  
  height = 10
)

library(RColorBrewer)

df_filtered = df_corrected[df_corrected$padj<0.05 & abs(df_corrected$log2FoldChange_std) > 0.58496250072, ]
region_counts <- df_allregions_RP_filtered %>%
  group_by(region) %>%
  summarise(
    total_rows = n(),
    unique_V1_count = n_distinct(V1)
  )

detection_counts <- df_filtered %>%
  group_by(region, V1) %>%
  summarise(
    detection_count = n(),
    avg_abs_log2FC = mean(abs(log2FoldChange_std))
  ) %>%
  ungroup()

gene_order <- detection_counts %>%
  group_by(V1) %>%
  summarise(total_detections = sum(detection_count)) %>%
  arrange(desc(total_detections)) %>%
  pull(V1) 
complete_counts <- detection_counts %>%
  complete(V1 = gene_order, region, fill = list(detection_count = 0))

complete_counts$V1 <- factor(complete_counts$V1, levels = gene_order)
final_plot_fixed_axes <- ggplot(complete_counts, aes(x = V1, y = detection_count)) +
  geom_col(aes(fill = region), show.legend = FALSE) +
  
  facet_wrap(~ region, ncol = 2) + 
  
  labs(
    title = NULL,
    x = "RP",
    y = "Number of DE Events"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.background = element_rect(fill = "darkslategray"),
    strip.text = element_text(color = "white", face = "bold")
  )

print(final_plot_fixed_axes)

dir.create("Figuras/Figure_S5", recursive = TRUE, showWarnings = FALSE)
ggsave("Figuras/Figure_S5/FigureS5_final_plot_fixed_axes.svg", plot = final_plot_fixed_axes, width = 18, height = 9)
ggsave("Figuras/Figure_S5/FigureS5_final_plot_fixed_axes.png", plot = final_plot_fixed_axes, width = 18, height = 9, dpi = 300)

