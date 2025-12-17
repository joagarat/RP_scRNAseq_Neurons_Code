
library(Matrix)
library(data.table)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(lme4)
library(lmerTest)
library(nlme)
library(readxl)
library(stringr)
library(ggpubr)
library(purrr)
library(grid)
library(svglite)
library(Seurat)
library(scCustomize)
library(RColorBrewer)
library(DESeq2)
library(Matrix.utils)
library(tibble)
library(BiocParallel)

#Rps27
Rps27_Gabavs_Glut <- read_excel(path = "data/Immuno_Quantification/Cuantificacion_Rps27_NeuN_Gad67.xlsx")  
Rps27_Gabavs_Glut = as.data.frame(Rps27_Gabavs_Glut)
Rps27_Gabavs_Glut <- Rps27_Gabavs_Glut %>%
  mutate(
    curated_cell_type = case_when(
      str_to_lower(cell_type) == "glut" ~ "Glut",
      str_to_lower(cell_type) == "gaba" ~ "GABA",
      TRUE ~ NA_character_
    )
  )
Rps27_Gabavs_Glut <- Rps27_Gabavs_Glut %>%
  filter(!is.na(curated_cell_type))
colnames(Rps27_Gabavs_Glut)[which(colnames(Rps27_Gabavs_Glut) == "Vol (pix)")] <- "Vol_pix"
colnames(Rps27_Gabavs_Glut)[which(colnames(Rps27_Gabavs_Glut) == "CZ (pix)")] <- "CZ_pix"
colnames(Rps27_Gabavs_Glut)[which(colnames(Rps27_Gabavs_Glut) == "S27 IntDen")] <- "S27_IntDen"
colnames(Rps27_Gabavs_Glut)[which(colnames(Rps27_Gabavs_Glut) == "S27 Mean")] <- "S27_Mean"
colnames(Rps27_Gabavs_Glut)

detectar_rango_z_valido <- function(df_imagen, 
                                    columna_intensidad = "S27_Mean", 
                                    columna_z = "CZ_pix", 
                                    threshold_ratio = 0.7) {
  df_binned <- df_imagen %>%
    group_by(.data[[columna_z]]) %>%
    summarise(mean_Mean = mean(.data[[columna_intensidad]], na.rm = TRUE), .groups = "drop")
  p95 <- quantile(df_binned$mean_Mean, probs = 0.95, na.rm = TRUE)
  top_mean <- df_binned %>%
    filter(mean_Mean >= p95) %>%
    summarise(mean_top = mean(mean_Mean)) %>%
    pull(mean_top)
  threshold_value <- threshold_ratio * top_mean
  z_range <- df_binned %>%
    filter(mean_Mean > threshold_value) %>%
    summarise(min_z = min(.data[[columna_z]]),
              max_z = max(.data[[columna_z]]))
  return(z_range)
}

rango_z_por_imagen <- Rps27_Gabavs_Glut %>%
  group_split(Img) %>%
  map_dfr(function(df_Img) {
    z_rango <- detectar_rango_z_valido(df_Img, threshold_ratio = 0.7) 
    tibble(
      Img = unique(df_Img$Img),
      min_z = z_rango$min_z,
      max_z = z_rango$max_z
    )
  })

Rps27_anotado <- Rps27_Gabavs_Glut %>%
  left_join(rango_z_por_imagen, by = "Img") %>%
  mutate(
    pasa_filtro_z = CZ_pix >= min_z & CZ_pix <= max_z
  )

Rps27_filtrado <- Rps27_anotado %>%
  filter(pasa_filtro_z)

imagenes = unique(Rps27_filtrado$Img)


df_plot <- Rps27_anotado %>% filter(Img == 1554)

p <- ggplot(df_plot, aes(x = CZ_pix, y = S27_Mean, color = pasa_filtro_z)) +
  geom_jitter(width = 0.3, alpha = 0.6) +
  geom_vline(xintercept = unique(df_plot$min_z), linetype = "dashed", color = "red") +
  geom_vline(xintercept = unique(df_plot$max_z), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "dodgerblue")) +
  labs(
    title = NULL,
    x = "Z plane (CZ_pix)",
    y = "Mean intensity (S27_Mean)",
    color = "Passes Z filter"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )
dir.create("Figuras/Figura_S7", recursive = TRUE, showWarnings = FALSE)
svglite(
  "Figuras/Figura_S7/Fig_S7A.svg",
  width = 2.5,
  height = 2.5
)

print(p)
dev.off()

Rps27_filtrado <- Rps27_filtrado %>%
  mutate(Img_type = interaction(Img, curated_cell_type),
         Img = as.factor(Img)) %>%
  arrange(Donor, Img, curated_cell_type) %>%
  mutate(Img_type = factor(Img_type, levels = unique(Img_type)))

datos_filtrados <- subset(Rps27_filtrado, Vol_pix >= 10000)

summary_df_filtrado <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  group_by(Donor, Img, curated_cell_type) %>%
  summarise(cell_count = n(), .groups = "drop")

percentages_df_filtrado <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  select(Img, Donor, curated_cell_type) %>%
  group_by(Donor, Img, curated_cell_type) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Donor, Img) %>%
  mutate(
    total_cells = sum(cell_count),
    percent = 100 * cell_count / total_cells
  ) %>%
  ungroup()


percentages_by_mouse <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  dplyr::count(Donor, curated_cell_type, name = "cell_count") %>% 
  group_by(Donor) %>%
  mutate(
    total_cells = sum(cell_count),
    raw_percent = 100 * cell_count / total_cells,   
    percent_2d  = round(raw_percent, 2)      
  ) %>%
  ungroup() %>%
  arrange(Donor, desc(cell_count))


library(ggplot2)

panel_percentage_plot <- ggplot(percentages_by_mouse, aes(x = Donor, y = raw_percent, fill = curated_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = NULL,
    x = "Mouse",
    y = "Neuron (%)",
    fill = "Cell Type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))   # legend in 2 rows

print(panel_percentage_plot)

ggsave(
  "Figuras/Figura_S7/Fig_S7B.svg",
  plot   = panel_percentage_plot,
  width  = 1.7,
  height = 2.5,
  units  = "in"
)

datos_filtrados <- datos_filtrados %>%
  group_by(Img) %>%
  mutate(
    IntDen_Norm = S27_IntDen / mean(S27_IntDen, na.rm = TRUE),
    Mean_Norm = S27_Mean / mean(S27_Mean, na.rm = TRUE)
  ) %>%
  ungroup()


datos_long_filtrados <- datos_filtrados %>%
  select(Donor, curated_cell_type,
         Norm_IntDen = `IntDen_Norm`,
         Norm_Mean   = `Mean_Norm`) %>%
  pivot_longer(cols = starts_with("Norm_"), names_to = "Medida", values_to = "Valor") %>%
  mutate(Medida = recode(Medida,
                         "Norm_IntDen" = "IntDen normalizado",
                         "Norm_Mean" = "Mean normalizado"))

modeloMean <- lme(S27_Mean ~ curated_cell_type, random = ~1 | Donor, data = datos_filtrados)
resumen <- capture.output(summary(modeloMean))

summary(modeloMean)

medias_por_raton <- datos_filtrados %>%
  group_by(curated_cell_type, Donor) %>%
  summarise(media = mean(Mean_Norm), .groups = "drop")
max_y <- max(datos_filtrados$Mean_Norm, na.rm = TRUE)
max_y <- max(datos_filtrados$Mean_Norm, na.rm = TRUE)

pos_barra <- max_y * 1.18        
pos_asterisco <- max_y * 1.20    
pos_texto_cambio <- max_y * 1.35 
mean_gaba <- mean(medias_por_raton$media[medias_por_raton$curated_cell_type == "GABA"], na.rm = TRUE)
mean_glut <- mean(medias_por_raton$media[medias_por_raton$curated_cell_type == "Glut"], na.rm = TRUE)
percent_change <- round(((mean_gaba - mean_glut) / mean_glut) * 100, 0)
panel_violin_plot_s27_protein <- ggplot(datos_filtrados, aes(x = curated_cell_type, y = Mean_Norm)) +
  geom_jitter(aes(color = curated_cell_type), width = 0.2, size = 0.5, alpha = 0.4) + 
  geom_violin(aes(color = curated_cell_type), fill = NA, size = 0.5, trim = FALSE) + 
  stat_summary(data = medias_por_raton, aes(y = media),
               fun = mean, geom = "crossbar", width = 0.3, 
               size = 0.2, color = "black") +
  geom_jitter(data = medias_por_raton, aes(y = media),
              width = 0.08, size = 1.5, color = "black", shape = 19) + 
  geom_segment(aes(x = 1, xend = 2, y = pos_barra, yend = pos_barra), color = "black", size = 0.3) + 
  geom_segment(aes(x = 1, xend = 1, y = pos_barra - (max_y*0.03), yend = pos_barra), color = "black", size = 0.3) + 
  geom_segment(aes(x = 2, xend = 2, y = pos_barra - (max_y*0.03), yend = pos_barra), color = "black", size = 0.3) + 
  annotate("text", x = 1.5, y = pos_asterisco, label = "***", size = 3, color = "black", vjust = 0) + 
  annotate("text", 
           x = 1.5, 
           y = pos_texto_cambio, 
           label = paste0(percent_change, "% Change"), 
           hjust = 0.5, 
           vjust = 1, 
           size = 2.5, fontface = "bold") + 
  labs(
    title = "Rps27 Protein levels", 
    y = "Normalized Protein Levels",
    x = NULL 
  ) +
  scale_color_manual(values = c("GABA" = "#F8766D", "Glut" = "#00BFC4")) + 
  scale_x_discrete(labels = c("GABA" = "GABAergic", "Glut" = "Glutamatergic")) +
  coord_cartesian(ylim = c(NA, 3.25)) +
  theme_classic() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), 
    axis.text = element_text(size = 8), 
    axis.title.y = element_text(size = 9), 
    axis.line = element_line(size = 0.3), 
    axis.ticks = element_line(size = 0.3)
  )

print(panel_violin_plot_s27_protein)
dir.create("Figuras/Figura_9", recursive = TRUE, showWarnings = FALSE)
ggsave(
  "Figuras/Figura_9/Figura_9C.png",
  plot = panel_violin_plot_s27_protein,
  width = 2.5, 
  height = 2.5, 
  units = "in",
  dpi = 300
)

load("data/Smartseq2dataset_seurat_filtered.RData")
DefaultAssay(seurat_merged) <- "RNA"
RP = read.csv("data/RP.csv")
RP = as.character(RP$x)
RP =gsub("Rack1", "Gnb2l1", x = RP)
gen_a_plotear <- "Rps27" 
datos_scrna <- FetchData(seurat_merged, 
                         vars = c(gen_a_plotear, "class_label", "donor_label"))
colnames(datos_scrna)[1] <- "Expression"
datos_scrna_filtrados <- datos_scrna %>%
  filter(class_label %in% c("GABAergic", "Glutamatergic"))
medias_por_donor_scrna <- datos_scrna_filtrados %>%
  group_by(donor_label, class_label) %>%
  summarise(media_expresion = mean(Expression, na.rm = TRUE)) %>%
  ungroup()
library(ggplot2)
altura_barra_scrna <- 2.2 
gen_de_interes <- "Rps27" 
counts <- GetAssayData(seurat_merged, assay = "RNA", slot = "counts")
metadata <- seurat_merged@meta.data %>%
  select(class_label, donor_label)
metadata$sample_id <- paste(metadata$class_label, metadata$donor_label, sep = "_")
pseudobulk_counts <- Matrix.utils::aggregate.Matrix(t(counts), 
                                                    groupings = metadata$sample_id, 
                                                    fun = "sum")
pseudobulk_counts <- t(pseudobulk_counts)
length(colnames(pseudobulk_counts))
coldata <- data.frame(sample = colnames(pseudobulk_counts)) %>%
  mutate(class_label = stringr::str_extract(sample, "GABAergic|Glutamatergic"),
         donor = stringr::str_extract(sample, "[^_]+$")) %>%
  column_to_rownames("sample")
coldata <- coldata[colnames(pseudobulk_counts), ]


dds <- DESeqDataSetFromMatrix(countData = pseudobulk_counts,
                              colData = coldata,
                              design = ~ class_label) 

register(MulticoreParam(10))
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, contrast = c("class_label", "Glutamatergic", "GABAergic"))
resultado_gen <- as.data.frame(res)[gen_de_interes, ]
mi_l2fc <- resultado_gen$log2FoldChange
mi_padj <- resultado_gen$padj
fc = 2^mi_l2fc
fc = 1/fc
datos_plot <- datos_scrna_filtrados %>% 
  filter(class_label %in% c("GABAergic", "Glutamatergic"))
medias_plot <- medias_por_donor_scrna %>% 
  filter(class_label %in% c("GABAergic", "Glutamatergic"))
max_y_mrna <- max(datos_plot$Expression, na.rm = TRUE)
pos_barra <- max_y_mrna * 1.12        
pos_asterisco <- max_y_mrna * 1.14    
pos_texto_cambio <- max_y_mrna * 1.35 
mean_gaba <- mean(medias_plot$media_expresion[medias_plot$class_label == "GABAergic"], na.rm = TRUE)
mean_glut <- mean(medias_plot$media_expresion[medias_plot$class_label == "Glutamatergic"], na.rm = TRUE)
percent_change <- round(100*(fc-1))
my_violin_plot_mrna <- ggplot(datos_plot, aes(x = class_label, y = Expression)) +
  geom_jitter(aes(color = class_label), width = 0.2, size = 0.5, alpha = 0.2) + 
  geom_violin(aes(color = class_label), fill = NA, size = 0.5, trim = FALSE) + 
  stat_summary(data = medias_plot, aes(y = media_expresion), fun = mean, 
               geom = "crossbar", width = 0.3, size = 0.2, color = "black") +
  geom_jitter(data = medias_plot, aes(y = media_expresion),
              width = 0.08, size = 1.5, color = "black", shape = 19) + 
  geom_segment(aes(x = 1, xend = 2, y = pos_barra, yend = pos_barra), color = "black", size = 0.3) +
  geom_segment(aes(x = 1, xend = 1, y = pos_barra - (max_y_mrna*0.03), yend = pos_barra), color = "black", size = 0.3) +
  geom_segment(aes(x = 2, xend = 2, y = pos_barra - (max_y_mrna*0.03), yend = pos_barra), color = "black", size = 0.3) +
  annotate("text", x = 1.5, y = pos_asterisco, label = "****", size = 3, color = "black", vjust = 0) +
  annotate("text", 
           x = 1.5, 
           y = pos_texto_cambio, 
           label = paste0(percent_change, "% Change"), 
           hjust = 0.5, 
           vjust = 1, 
           size = 2.5, fontface = "bold") +
  labs(
    title = "Rps27 mRNA levels",
    y = "Normalized mRNA Levels",
    x = NULL 
  ) +
  scale_color_manual(values = c("GABAergic" = "#F8766D", "Glutamatergic" = "#00BFC4")) +
  coord_cartesian(ylim = c(NA, max_y_mrna * 1.45)) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3)
  )

print(my_violin_plot_mrna)

ggsave(
  "Figuras/Figura_9/Figura_9B.png",
  plot = my_violin_plot_mrna,
  width = 2.5, 
  height = 2.5, 
  units = "in",
  dpi = 300
)

#Rpl29
Rpl29_Gabavs_Glut <- read_excel(path = "data/Immuno_Quantification/Cuantificacion_Rpl29_NeuN_Gad67.xlsx")  
Rpl29_Gabavs_Glut = as.data.frame(Rpl29_Gabavs_Glut)
Rpl29_Gabavs_Glut <- Rpl29_Gabavs_Glut %>%
  mutate(
    curated_cell_type = case_when(
      str_to_lower(cell_type) == "glut" ~ "Glut",
      str_to_lower(cell_type) == "gaba" ~ "GABA",
      TRUE ~ NA_character_
    )
  )
Rpl29_Gabavs_Glut <- Rpl29_Gabavs_Glut %>%
  filter(!is.na(curated_cell_type))
colnames(Rpl29_Gabavs_Glut)[which(colnames(Rpl29_Gabavs_Glut) == "Vol (pix)")] <- "Vol_pix"
colnames(Rpl29_Gabavs_Glut)[which(colnames(Rpl29_Gabavs_Glut) == "CZ (pix)")] <- "CZ_pix"
colnames(Rpl29_Gabavs_Glut)[which(colnames(Rpl29_Gabavs_Glut) == "L29 IntDen")] <- "L29_IntDen"
colnames(Rpl29_Gabavs_Glut)[which(colnames(Rpl29_Gabavs_Glut) == "L29 Mean")] <- "L29_Mean"
colnames(Rpl29_Gabavs_Glut)
detectar_rango_z_valido <- function(df_imagen, 
                                    columna_intensidad = "L29_Mean", 
                                    columna_z = "CZ_pix", 
                                    threshold_ratio = 0.7) {
  df_binned <- df_imagen %>%
    group_by(.data[[columna_z]]) %>%
    summarise(mean_Mean = mean(.data[[columna_intensidad]], na.rm = TRUE), .groups = "drop")
  p95 <- quantile(df_binned$mean_Mean, probs = 0.95, na.rm = TRUE)
  top_mean <- df_binned %>%
    filter(mean_Mean >= p95) %>%
    summarise(mean_top = mean(mean_Mean)) %>%
    pull(mean_top)
  threshold_value <- threshold_ratio * top_mean
  z_range <- df_binned %>%
    filter(mean_Mean > threshold_value) %>%
    summarise(min_z = min(.data[[columna_z]]),
              max_z = max(.data[[columna_z]]))
  return(z_range)
}

rango_z_por_imagen <- Rpl29_Gabavs_Glut %>%
  group_split(img) %>%
  map_dfr(function(df_img) {
    z_rango <- detectar_rango_z_valido(df_img, threshold_ratio = 0.7) 
    tibble(
      img = unique(df_img$img),
      min_z = z_rango$min_z,
      max_z = z_rango$max_z
    )
  })

Rpl29_anotado <- Rpl29_Gabavs_Glut %>%
  left_join(rango_z_por_imagen, by = "img") %>%
  mutate(
    pasa_filtro_z = CZ_pix >= min_z & CZ_pix <= max_z
  )

Rpl29_filtrado <- Rpl29_anotado %>%
  filter(pasa_filtro_z)

summary_df <- Rpl29_filtrado %>%
  filter(!is.na(curated_cell_type)) %>%
  group_by(Mouse, img, curated_cell_type) %>%
  summarise(cell_count = n(), .groups = "drop")

Rpl29_filtrado <- Rpl29_filtrado %>%
  mutate(img_type = interaction(img, curated_cell_type),
         img = as.factor(img)) %>%
  arrange(Mouse, img, curated_cell_type) %>%
  mutate(img_type = factor(img_type, levels = unique(img_type)))

Rpl29_filtrado <- Rpl29_filtrado %>%
  group_by(img) %>%
  mutate(
    IntDen_Norm = L29_IntDen / mean(L29_IntDen, na.rm = TRUE),
    Mean_Norm = L29_Mean / mean(L29_Mean, na.rm = TRUE)
  ) %>%
  ungroup()

library(gridExtra)

umbrales <- c(5000, 10000, 15000, 20000)
datos_filtrados <- subset(Rpl29_filtrado, Vol_pix >= 10000)

summary_df_filtrado <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  group_by(Mouse, img, curated_cell_type) %>%
  summarise(cell_count = n(), .groups = "drop")

percentages_df_filtrado <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  select(img, Mouse, curated_cell_type) %>%
  group_by(Mouse, img, curated_cell_type) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Mouse, img) %>%
  mutate(
    total_cells = sum(cell_count),
    percent = 100 * cell_count / total_cells
  ) %>%
  ungroup()

library(ggplot2)
library(svglite)

percentages_by_mouse <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  dplyr::count(Mouse, curated_cell_type, name = "cell_count") %>% 
  group_by(Mouse) %>%
  mutate(
    total_cells = sum(cell_count),
    raw_percent = 100 * cell_count / total_cells,    
    percent_2d  = round(raw_percent, 2)            
  ) %>%
  ungroup() %>%
  arrange(Mouse, desc(cell_count))

panel_percentage_plot <- ggplot(percentages_by_mouse, aes(x = Mouse, y = raw_percent, fill = curated_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = NULL,
    x = "Mouse",
    y = "Neuron %",
    fill = "Cell Type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))  
print(panel_percentage_plot)

ggsave(
  "Figuras/Figura_S7/Fig_S7C.svg",
  plot   = panel_percentage_plot,
  width  = 1.7,
  height = 2.5,
  units  = "in"
)

datos_filtrados <- datos_filtrados %>%
  group_by(img) %>%
  mutate(
    IntDen_Norm = L29_IntDen / mean(L29_IntDen, na.rm = TRUE),
    Mean_Norm = L29_Mean / mean(L29_Mean, na.rm = TRUE)
  ) %>%
  ungroup()

datos_long_filtrados <- datos_filtrados %>%
  select(Mouse, curated_cell_type,
         Norm_IntDen = `IntDen_Norm`,
         Norm_Mean   = `Mean_Norm`) %>%
  pivot_longer(cols = starts_with("Norm_"), names_to = "Medida", values_to = "Valor") %>%
  mutate(Medida = recode(Medida,
                         "Norm_IntDen" = "IntDen normalizado",
                         "Norm_Mean" = "Mean normalizado"))


modeloMean <- lme(L29_Mean ~ curated_cell_type, random = ~1 | Mouse, data = datos_filtrados)
resumen <- capture.output(summary(modeloMean))

medias_por_raton <- datos_filtrados %>%
  group_by(curated_cell_type, Mouse) %>%
  summarise(media = mean(Mean_Norm), .groups = "drop")

panel_violin_plot_rpl29 <- ggplot(datos_filtrados, aes(x = curated_cell_type, y = Mean_Norm)) +
  geom_jitter(aes(color = curated_cell_type), width = 0.2, size = 0.5, alpha = 0.4) +
  geom_violin(aes(color = curated_cell_type), fill = NA, size = 0.5, trim = FALSE) +
  stat_summary(data = medias_por_raton, aes(y = media),
               fun = mean, geom = "crossbar", width = 0.3,
               size = 0.2, color = "black") +
  geom_jitter(data = medias_por_raton, aes(y = media),
              width = 0.08, size = 1.5, color = "black", shape = 19) +
  labs(
    title = "Rpl29 Protein levels",
    y = "Normalized Protein Levels",
    x = NULL
  ) +
  scale_color_manual(values = c("GABA" = "#F8766D", "Glut" = "#00BFC4")) +
  scale_x_discrete(labels = c("GABA" = "GABAergic", "Glut" = "Glutamatergic")) +
  coord_cartesian(ylim = c(NA, 3.25)) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3)
  )

print(panel_violin_plot_rpl29)

ggsave(
  "Figuras/Figura_S7/Fig_S7H.svg",
  plot   = panel_violin_plot_rpl29,
  width  = 2.5,
  height = 2.5,
  units  = "in",
  dpi    = 300
)

gen_a_plotear <- "Rpl29" 
datos_scrna <- FetchData(seurat_merged, 
                         vars = c(gen_a_plotear, "class_label", "donor_label"))
colnames(datos_scrna)[1] <- "Expression"
datos_scrna_filtrados <- datos_scrna %>%
  filter(class_label %in% c("GABAergic", "Glutamatergic"))
medias_por_donor_scrna <- datos_scrna_filtrados %>%
  group_by(donor_label, class_label) %>%
  summarise(media_expresion = mean(Expression, na.rm = TRUE)) %>%
  ungroup()

library(ggplot2)

altura_barra_scrna <- 2.2

library(ggplot2)
library(dplyr)

gen_de_interes <- "Rpl29" 
counts <- GetAssayData(seurat_merged, assay = "RNA", slot = "counts")
metadata <- seurat_merged@meta.data %>%
  select(class_label, donor_label)
metadata$sample_id <- paste(metadata$class_label, metadata$donor_label, sep = "_")
pseudobulk_counts <- Matrix.utils::aggregate.Matrix(t(counts), 
                                                    groupings = metadata$sample_id, 
                                                    fun = "sum")
pseudobulk_counts <- t(pseudobulk_counts)
length(colnames(pseudobulk_counts))
coldata <- data.frame(sample = colnames(pseudobulk_counts)) %>%
  mutate(class_label = stringr::str_extract(sample, "GABAergic|Glutamatergic"),
         donor = stringr::str_extract(sample, "[^_]+$")) %>%
  column_to_rownames("sample")
coldata <- coldata[colnames(pseudobulk_counts), ]
dds <- DESeqDataSetFromMatrix(countData = pseudobulk_counts,
                              colData = coldata,
                              design = ~ class_label)

register(MulticoreParam(10))
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, contrast = c("class_label", "Glutamatergic", "GABAergic"))
resultado_gen <- as.data.frame(res)[gen_de_interes, ]
mi_l2fc <- resultado_gen$log2FoldChange
mi_padj <- resultado_gen$padj
fc = 2^mi_l2fc
fc = 1/fc
datos_plot <- datos_scrna_filtrados %>% 
  filter(class_label %in% c("GABAergic", "Glutamatergic"))
medias_plot <- medias_por_donor_scrna %>% 
  filter(class_label %in% c("GABAergic", "Glutamatergic"))

max_y_mrna <- max(datos_plot$Expression, na.rm = TRUE)
pos_barra <- max_y_mrna * 1.12        
pos_asterisco <- max_y_mrna * 1.14    
pos_texto_cambio <- max_y_mrna * 1.35 

mean_gaba <- mean(medias_plot$media_expresion[medias_plot$class_label == "GABAergic"], na.rm = TRUE)
mean_glut <- mean(medias_plot$media_expresion[medias_plot$class_label == "Glutamatergic"], na.rm = TRUE)
percent_change <- round(100*(fc-1))

my_violin_plot_mrna <- ggplot(datos_plot, aes(x = class_label, y = Expression)) +
  geom_jitter(aes(color = class_label), width = 0.2, size = 0.5, alpha = 0.2) + 
  geom_violin(aes(color = class_label), fill = NA, size = 0.5, trim = FALSE) + 
  stat_summary(data = medias_plot, aes(y = media_expresion), fun = mean, 
               geom = "crossbar", width = 0.3, size = 0.2, color = "black") +
  geom_jitter(data = medias_plot, aes(y = media_expresion),
              width = 0.08, size = 1.5, color = "black", shape = 19) + 
  geom_segment(aes(x = 1, xend = 2, y = pos_barra, yend = pos_barra), color = "black", size = 0.3) +
  geom_segment(aes(x = 1, xend = 1, y = pos_barra - (max_y_mrna*0.03), yend = pos_barra), color = "black", size = 0.3) +
  geom_segment(aes(x = 2, xend = 2, y = pos_barra - (max_y_mrna*0.03), yend = pos_barra), color = "black", size = 0.3) +
  annotate("text", x = 1.5, y = pos_asterisco, label = "****", size = 3, color = "black", vjust = 0) +
  annotate("text", 
           x = 1.5, 
           y = pos_texto_cambio, 
           label = paste0(percent_change, "% Change"), 
           hjust = 0.5, 
           vjust = 1, 
           size = 2.5, fontface = "bold") +
  labs(
    title = "Rpl29 mRNA levels", # Título ajustado
    y = "Normalized mRNA Levels",
    x = NULL 
  ) +
  scale_color_manual(values = c("GABAergic" = "#F8766D", "Glutamatergic" = "#00BFC4")) +
  coord_cartesian(ylim = c(NA, max_y_mrna * 1.45)) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3)
  )

print(my_violin_plot_mrna)

ggsave(
  "Figuras/Figura_S7/Figura_S7G.png",
  plot = my_violin_plot_mrna,
  width = 2.5, 
  height = 2.5, 
  units = "in",
  dpi = 300
)

#Rpl21
Rpl21_Gabavs_Glut <- read_excel(path = "data/Immuno_Quantification/Cuantificacion_Rpl21_NeuN_Gad67.xlsx")  
Rpl21_Gabavs_Glut = as.data.frame(Rpl21_Gabavs_Glut)
Rpl21_Gabavs_Glut <- Rpl21_Gabavs_Glut %>%
  mutate(
    curated_cell_type = case_when(
      str_to_lower(cell_type) == "glut" ~ "Glut",
      str_to_lower(cell_type) == "gaba" ~ "GABA",
      TRUE ~ NA_character_
    )
  )
Rpl21_Gabavs_Glut <- Rpl21_Gabavs_Glut %>%
  filter(!is.na(curated_cell_type))
colnames(Rpl21_Gabavs_Glut)[which(colnames(Rpl21_Gabavs_Glut) == "Vol (pix)")] <- "Vol_pix"
colnames(Rpl21_Gabavs_Glut)[which(colnames(Rpl21_Gabavs_Glut) == "CZ (pix)")] <- "CZ_pix"
colnames(Rpl21_Gabavs_Glut)[which(colnames(Rpl21_Gabavs_Glut) == "L21 IntDen")] <- "L21_IntDen"
colnames(Rpl21_Gabavs_Glut)[which(colnames(Rpl21_Gabavs_Glut) == "L21 Mean")] <- "L21_Mean"
colnames(Rpl21_Gabavs_Glut)
detectar_rango_z_valido <- function(df_imagen, 
                                    columna_intensidad = "L21_Mean", 
                                    columna_z = "CZ_pix", 
                                    threshold_ratio = 0.7) {
  df_binned <- df_imagen %>%
    group_by(.data[[columna_z]]) %>%
    summarise(mean_Mean = mean(.data[[columna_intensidad]], na.rm = TRUE), .groups = "drop")
  p95 <- quantile(df_binned$mean_Mean, probs = 0.95, na.rm = TRUE)
  top_mean <- df_binned %>%
    filter(mean_Mean >= p95) %>%
    summarise(mean_top = mean(mean_Mean)) %>%
    pull(mean_top)
  threshold_value <- threshold_ratio * top_mean
  z_range <- df_binned %>%
    filter(mean_Mean > threshold_value) %>%
    summarise(min_z = min(.data[[columna_z]]),
              max_z = max(.data[[columna_z]]))
  return(z_range)
}

rango_z_por_imagen <- Rpl21_Gabavs_Glut %>%
  group_split(img) %>%
  map_dfr(function(df_img) {
    z_rango <- detectar_rango_z_valido(df_img, threshold_ratio = 0.7) 
    tibble(
      img = unique(df_img$img),
      min_z = z_rango$min_z,
      max_z = z_rango$max_z
    )
  })

Rpl21_anotado <- Rpl21_Gabavs_Glut %>%
  left_join(rango_z_por_imagen, by = "img") %>%
  mutate(
    pasa_filtro_z = CZ_pix >= min_z & CZ_pix <= max_z
  )

Rpl21_filtrado <- Rpl21_anotado %>%
  filter(pasa_filtro_z)
imagenes = unique(Rpl21_filtrado$img)

Rpl21_filtrado <- Rpl21_filtrado %>%
  mutate(img_type = interaction(img, curated_cell_type),
         img = as.factor(img)) %>%
  arrange(Mouse, img, curated_cell_type) %>%
  mutate(img_type = factor(img_type, levels = unique(img_type)))

Rpl21_filtrado <- Rpl21_filtrado %>%
  group_by(img) %>%
  mutate(
    IntDen_Norm = L21_IntDen / mean(L21_IntDen, na.rm = TRUE),
    Mean_Norm = L21_Mean / mean(L21_Mean, na.rm = TRUE)
  ) %>%
  ungroup()

datos_long <- Rpl21_filtrado %>%
  select(Mouse, curated_cell_type,
         Norm_IntDen = `IntDen_Norm`,
         Norm_Mean   = `Mean_Norm`) %>%
  pivot_longer(cols = starts_with("Norm_"), names_to = "Medida", values_to = "Valor") %>%
  mutate(Medida = recode(Medida,
                         "Norm_IntDen" = "IntDen normalizado",
                         "Norm_Mean" = "Mean normalizado"))

umbrales <- c(5000, 10000, 15000, 20000)
datos_filtrados <- subset(Rpl21_filtrado, Vol_pix >= 10000)
summary_df_filtrado <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  group_by(Mouse, img, curated_cell_type) %>%
  summarise(cell_count = n(), .groups = "drop")
percentages_df_filtrado <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  select(img, Mouse, curated_cell_type) %>%
  group_by(Mouse, img, curated_cell_type) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Mouse, img) %>%
  mutate(
    total_cells = sum(cell_count),
    percent = 100 * cell_count / total_cells
  ) %>%
  ungroup()
percentages_by_mouse <- datos_filtrados %>%
  filter(!is.na(curated_cell_type)) %>%
  dplyr::count(Mouse, curated_cell_type, name = "cell_count") %>%  
  group_by(Mouse) %>%
  mutate(
    total_cells = sum(cell_count),
    raw_percent = 100 * cell_count / total_cells,    
    percent_2d  = round(raw_percent, 2)           
  ) %>%
  ungroup() %>%
  arrange(Mouse, desc(cell_count))

panel_percentage_by_mouse <- ggplot(percentages_by_mouse, aes(x = Mouse, y = raw_percent, fill = curated_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = NULL,
    x = "Mouse",
    y = "Neuron %",
    fill = "Cell Type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) 
print(panel_percentage_by_mouse)

ggsave(
  "Figuras/Figura_S7/Fig_S7D.svg",
  plot   = panel_percentage_by_mouse,
  width  = 1.7,
  height = 2.5,
  units  = "in"
)

datos_filtrados <- datos_filtrados %>%
  group_by(img) %>%
  mutate(
    IntDen_Norm = L21_IntDen / mean(L21_IntDen, na.rm = TRUE),
    Mean_Norm = L21_Mean / mean(L21_Mean, na.rm = TRUE)
  ) %>%
  ungroup()

datos_long_filtrados <- datos_filtrados %>%
  select(Mouse, curated_cell_type,
         Norm_IntDen = `IntDen_Norm`,
         Norm_Mean   = `Mean_Norm`) %>%
  pivot_longer(cols = starts_with("Norm_"), names_to = "Medida", values_to = "Valor") %>%
  mutate(Medida = recode(Medida,
                         "Norm_IntDen" = "IntDen normalizado",
                         "Norm_Mean" = "Mean normalizado"))

modeloMean <- lme(L21_Mean ~ curated_cell_type, random = ~1 | Mouse, data = datos_filtrados)
resumen <- capture.output(summary(modeloMean))
summary(modeloMean)

medias_por_raton <- datos_filtrados %>%
  group_by(curated_cell_type, Mouse) %>%
  summarise(media = mean(Mean_Norm), .groups = "drop")

panel_violin_plot_rpl21 <- ggplot(datos_filtrados, aes(x = curated_cell_type, y = Mean_Norm)) +
  geom_jitter(aes(color = curated_cell_type), width = 0.2, size = 0.5, alpha = 0.4) +
  geom_violin(aes(color = curated_cell_type), fill = NA, size = 0.5, trim = FALSE) +
  stat_summary(data = medias_por_raton, aes(y = media),
               fun = mean, geom = "crossbar", width = 0.3,
               size = 0.2, color = "black") +
  geom_jitter(data = medias_por_raton, aes(y = media),
              width = 0.08, size = 1.5, color = "black", shape = 19) +
  labs(
    title = "Rpl21 Protein levels",
    y = "Normalized Protein Levels",
    x = NULL
  ) +
  scale_color_manual(values = c("GABA" = "#F8766D", "Glut" = "#00BFC4")) +
  scale_x_discrete(labels = c("GABA" = "GABAergic", "Glut" = "Glutamatergic")) +
  coord_cartesian(ylim = c(NA, 2.5)) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3)
  )

print(panel_violin_plot_rpl21)

ggsave(
  "Figuras/Figura_S7/Fig_S7F.svg",
  plot   = panel_violin_plot_rpl21,
  width  = 2.5,
  height = 2.5,
  units  = "in",
  dpi    = 300
)

gen_a_plotear <- "Rpl21" 

datos_scrna <- FetchData(seurat_merged, 
                         vars = c(gen_a_plotear, "class_label", "donor_label"))

colnames(datos_scrna)[1] <- "Expression"
datos_scrna_filtrados <- datos_scrna %>%
  filter(class_label %in% c("GABAergic", "Glutamatergic"))
medias_por_donor_scrna <- datos_scrna_filtrados %>%
  group_by(donor_label, class_label) %>%
  summarise(media_expresion = mean(Expression, na.rm = TRUE)) %>%
  ungroup()
altura_barra_scrna <- 2.2
gen_de_interes <- "Rpl21" 
counts <- GetAssayData(seurat_merged, assay = "RNA", slot = "counts")
metadata <- seurat_merged@meta.data %>%
  select(class_label, donor_label)
metadata$sample_id <- paste(metadata$class_label, metadata$donor_label, sep = "_")
pseudobulk_counts <- Matrix.utils::aggregate.Matrix(t(counts), 
                                                    groupings = metadata$sample_id, 
                                                    fun = "sum")
pseudobulk_counts <- t(pseudobulk_counts)
length(colnames(pseudobulk_counts))
coldata <- data.frame(sample = colnames(pseudobulk_counts)) %>%
  mutate(class_label = stringr::str_extract(sample, "GABAergic|Glutamatergic"),
         donor = stringr::str_extract(sample, "[^_]+$")) %>%
  column_to_rownames("sample")
coldata <- coldata[colnames(pseudobulk_counts), ]

dds <- DESeqDataSetFromMatrix(countData = pseudobulk_counts,
                              colData = coldata,
                              design = ~ class_label) 

register(MulticoreParam(10))
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, contrast = c("class_label", "Glutamatergic", "GABAergic"))
resultado_gen <- as.data.frame(res)[gen_de_interes, ]
mi_l2fc <- resultado_gen$log2FoldChange
mi_padj <- resultado_gen$padj
fc = 2^mi_l2fc
fc = 1/fc
datos_plot <- datos_scrna_filtrados %>% 
  filter(class_label %in% c("GABAergic", "Glutamatergic"))

medias_plot <- medias_por_donor_scrna %>% 
  filter(class_label %in% c("GABAergic", "Glutamatergic"))

max_y_mrna <- max(datos_plot$Expression, na.rm = TRUE)

pos_barra <- max_y_mrna * 1.12        
pos_asterisco <- max_y_mrna * 1.14    
pos_texto_cambio <- max_y_mrna * 1.35 

mean_gaba <- mean(medias_plot$media_expresion[medias_plot$class_label == "GABAergic"], na.rm = TRUE)
mean_glut <- mean(medias_plot$media_expresion[medias_plot$class_label == "Glutamatergic"], na.rm = TRUE)
percent_change <- round(100*(fc-1))

my_violin_plot_mrna <- ggplot(datos_plot, aes(x = class_label, y = Expression)) +
  geom_jitter(aes(color = class_label), width = 0.2, size = 0.5, alpha = 0.2) + 
  geom_violin(aes(color = class_label), fill = NA, size = 0.5, trim = FALSE) + 
  stat_summary(data = medias_plot, aes(y = media_expresion), fun = mean, 
               geom = "crossbar", width = 0.3, size = 0.2, color = "black") +
  geom_jitter(data = medias_plot, aes(y = media_expresion),
              width = 0.08, size = 1.5, color = "black", shape = 19) + 
  geom_segment(aes(x = 1, xend = 2, y = pos_barra, yend = pos_barra), color = "black", size = 0.3) +
  geom_segment(aes(x = 1, xend = 1, y = pos_barra - (max_y_mrna*0.03), yend = pos_barra), color = "black", size = 0.3) +
  geom_segment(aes(x = 2, xend = 2, y = pos_barra - (max_y_mrna*0.03), yend = pos_barra), color = "black", size = 0.3) +
  annotate("text", x = 1.5, y = pos_asterisco, label = "****", size = 3, color = "black", vjust = 0) +
  annotate("text", 
           x = 1.5, 
           y = pos_texto_cambio, 
           label = paste0(percent_change, "% Change"), 
           hjust = 0.5, 
           vjust = 1, 
           size = 2.5, fontface = "bold") +
  labs(
    title = "Rpl21 mRNA levels", # Título ajustado
    y = "Normalized mRNA Levels",
    x = NULL 
  ) +
  scale_color_manual(values = c("GABAergic" = "#F8766D", "Glutamatergic" = "#00BFC4")) +
  coord_cartesian(ylim = c(NA, max_y_mrna * 1.45)) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3)
  )

print(my_violin_plot_mrna)

ggsave(
  "Figuras/Figura_S7/Figura_S7E.png",
  plot = my_violin_plot_mrna,
  width = 2.5, 
  height = 2.5, 
  units = "in",
  dpi = 300
)
