library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(ggrepel)
library(stringr)

RP <- read.csv("data/RP.csv")
RP <- as.character(RP$x)
RP <- gsub("Rack1", "Gnb2l1", x = RP)

inhibitory <- c("Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip")

load("data/Smartseq2dataset_seurat_filtered.RData")

seurat_split <- SplitObject(seurat_merged, split.by = "region_label")
seurat_split <- lapply(seurat_split, function(x) { DefaultAssay(x) <- "RNA"; x })

subclass_pass <- list(
  c("L5 PT CTX","L5 IT CTX","L4/5 IT CTX","L6 IT CTX","L6 CT CTX","L5 NP CTX","Pvalb","Vip","L2/3 IT CTX","Lamp5","Sst","Sst Chodl","Sncg","Car3","L6b CTX"),
  c("L4/5 IT CTX","L2/3 IT CTX","Lamp5","L5 IT CTX","L6 IT CTX","Vip","L6 CT CTX","L5 PT CTX","Car3","L5 NP CTX"),
  c("Car3","L2/3 IT CTX"),
  c("Vip","L5 NP CTX","L4/5 IT CTX","Sst","L6 IT CTX","L5 IT CTX","Lamp5","L2/3 IT CTX","L6 CT CTX","L6b CTX"),
  c("L5 PT CTX","L4 RSP-ACA"),
  c("L2/3 IT CTX","L6 IT CTX"),
  c("Sst","L5 NP CTX","L5 IT CTX","Vip","Lamp5","L4/5 IT CTX","Pvalb","L2/3 IT CTX","L6 IT CTX","L6 CT CTX","L6b CTX","Sncg"),
  c("DG","CA1-ProS","CA3","Lamp5","Vip","Sncg","Pvalb","Sst"),
  c("Sst","Sst Chodl","Sncg","Vip","L6 CT CTX","L4/5 IT CTX","L5 IT CTX","L2/3 IT CTX","Lamp5","L5 NP CTX","Pvalb","L6 IT CTX","L5 PT CTX"),
  c("L2/3 IT CTX","L4/5 IT CTX"),
  c("L4/5 IT CTX","L2/3 IT CTX","L5 IT CTX","L2/3 IT PPP","L6 CT CTX","L5 NP CTX","L6 IT CTX","L6b CTX","Sst","Lamp5","Vip","Pvalb")
)
regions <- c("VISp","VIS","SSs","SSp","RSPv","RSP","MOp","HIP","ALM","AI","ACA")

m <- matrix(nrow = 84, ncol = 88)
colnames(m) <- unlist(lapply(seq_along(subclass_pass), function(i) paste(subclass_pass[[i]], regions[i], sep = "_")))
rownames(m) <- RP

z <- 0
for (i in seq(1, 11)) {
  metadata_i <- as.data.frame(seurat_split[[i]]@meta.data)
  data       <- as.data.frame(GetAssayData(object = seurat_split[[i]], layer = "data"))
  for (p in seq_len(length(subclass_pass[[i]]))) {
    z <- z + 1
    cells       <- metadata_i[metadata_i$subclass_label %in% subclass_pass[[i]][p], ]
    datacellRP  <- data[RP, colnames(data) %in% rownames(cells)]
    m[, z]      <- rowMeans(datacellRP)
  }
}

remove(seurat_merged, seurat_split); gc()
RP <- gsub("Gnb2l1", "Rack1", x = RP)
subtypes_ss <- unique(gsub("_.*", "", colnames(m)))
m_ss <- sapply(subtypes_ss, function(st) {
  cols <- grep(paste0("^", st, "_"), colnames(m), value = TRUE)
  rowMeans(m[, cols, drop = FALSE], na.rm = TRUE)
})

seurat_adult <- readRDS("data/Normalized_Classified_Seurat_Adult")
seurat_aged  <- readRDS("data/Normalized_Classified_Seurat_Aged")

calculate_avg_expression_matrix <- function(seurat_object, rp_genes) {
  avg_expr <- AverageExpression(
    object   = seurat_object,
    features = rp_genes,
    group.by = "predicted.id",
    assay    = "RNA"
  )
  return(avg_expr$RNA)
}

m_adult <- calculate_avg_expression_matrix(seurat_adult, RP)
m_aged  <- calculate_avg_expression_matrix(seurat_aged,  RP)
remove(seurat_adult, seurat_aged); gc()

common_genes   <- intersect(rownames(m_adult), rownames(m_aged))
common_subtypes <- intersect(colnames(m_adult), colnames(m_aged))
m_adult <- m_adult[common_genes, common_subtypes]
m_aged  <- m_aged[common_genes,  common_subtypes]
colnames(m_adult) <- paste0(colnames(m_adult), "_adult")
colnames(m_aged)  <- paste0(colnames(m_aged),  "_aged")
m_combined_aging  <- cbind(m_adult, m_aged) 

seurat_stress <- LoadSeuratRds("data/Hing_2024_PFC_stress/SeuratProject.Rds")

seurat_stress$stress_group <- str_extract(seurat_stress$replica, "^[^_]+")
seurat_low    <- subset(seurat_stress, subset = stress_group == "low")
seurat_medium <- subset(seurat_stress, subset = stress_group == "medium")
seurat_high   <- subset(seurat_stress, subset = stress_group == "high")
remove(seurat_stress); gc()

m_low    <- calculate_avg_expression_matrix(seurat_low,    RP)
m_medium <- calculate_avg_expression_matrix(seurat_medium, RP)
m_high   <- calculate_avg_expression_matrix(seurat_high,   RP)
remove(seurat_low, seurat_medium, seurat_high); gc()

common_genes    <- intersect(rownames(m_low), intersect(rownames(m_medium), rownames(m_high)))
common_subtypes <- intersect(colnames(m_low), intersect(colnames(m_medium), colnames(m_high)))
m_low    <- m_low[common_genes,    common_subtypes]
m_medium <- m_medium[common_genes, common_subtypes]
m_high   <- m_high[common_genes,   common_subtypes]
colnames(m_low)    <- paste0(colnames(m_low),    "_low")
colnames(m_medium) <- paste0(colnames(m_medium), "_medium")
colnames(m_high)   <- paste0(colnames(m_high),   "_high")
m_combined_stress  <- cbind(m_low, m_medium, m_high)  

calc_zscore_diff <- function(avg_mat, dataset_label,
                             genes_highlight = c("Rpl21", "Rps27")) {
  
  mat_scaled <- scale(as.matrix(avg_mat))
  mat_scaled[is.na(mat_scaled)] <- 0
  
  colnames_clean <- gsub("_adult$|_aged$|_low$|_medium$|_high$", "", colnames(mat_scaled))
  gaba_cols <- colnames(mat_scaled)[colnames_clean %in% inhibitory]
  glut_cols <- colnames(mat_scaled)[!colnames_clean %in% inhibitory]
  
  cat(dataset_label, " GABA cols:", length(gaba_cols),
      "| Glut cols:", length(glut_cols), "\n")
  
  mean_gaba_z <- rowMeans(mat_scaled[, gaba_cols, drop = FALSE], na.rm = TRUE)
  mean_glut_z <- rowMeans(mat_scaled[, glut_cols, drop = FALSE], na.rm = TRUE)
  diff_z      <- mean_gaba_z - mean_glut_z
  
  pvals <- sapply(rownames(mat_scaled), function(gene) {
    gaba_vals <- mat_scaled[gene, gaba_cols]
    glut_vals <- mat_scaled[gene, glut_cols]
    if (length(gaba_vals) >= 2 && length(glut_vals) >= 2) {
      tryCatch(t.test(gaba_vals, glut_vals)$p.value, error = function(e) NA)
    } else NA
  })
  
  df <- data.frame(
    Gene         = rownames(mat_scaled),
    mean_GABA_z  = mean_gaba_z,
    mean_Glut_z  = mean_glut_z,
    diff_z       = diff_z,
    pval         = pvals,
    is_highlight = rownames(mat_scaled) %in% genes_highlight
  ) %>%
    arrange(desc(diff_z)) %>%
    mutate(
      rank    = row_number(),
      Dataset = dataset_label,
      padj    = p.adjust(pval, method = "BH"),
      sig     = ifelse(!is.na(padj) & padj < 0.05, "sig", "ns")
    )
  
  return(df)
}

diff_smartseq <- calc_zscore_diff(m_ss,             "Smart-seq")
diff_aging    <- calc_zscore_diff(m_combined_aging,  "Jin et al. - Aging")
diff_stress   <- calc_zscore_diff(m_combined_stress, "Hing et al. - Stress")

for (df in list(diff_smartseq, diff_aging, diff_stress)) {
  ds <- unique(df$Dataset)
  cat("\n===", ds, "===\n")
  cat("Top 10 genes con mayor z-score en GABAérgicas:\n")
  print(df %>% arrange(desc(diff_z)) %>%
          select(Gene, diff_z, padj, sig) %>% slice(1:10))
  for (g in c("Rpl21", "Rps27")) {
    row <- df[df$Gene == g, ]
    cat(g, " diff_z:", round(row$diff_z, 3),
        "| rank:", row$rank, "de", nrow(df),
        "| padj:", round(row$padj, 4),
        "| sig:", row$sig, "\n")
  }
}

library(ggrepel)

df_concordance <- diff_smartseq %>%
  select(Gene, diff_z_smartseq = diff_z, is_highlight) %>%
  left_join(
    diff_aging %>% select(Gene, diff_z_aging = diff_z),
    by = "Gene"
  ) %>%
  left_join(
    diff_stress %>% select(Gene, diff_z_stress = diff_z),
    by = "Gene"
  ) %>%
  mutate(
    concordant_aging = sign(diff_z_smartseq) == sign(diff_z_aging),
    concordant_stress = sign(diff_z_smartseq) == sign(diff_z_stress),
    label = ifelse(is_highlight, Gene, "")
  )

make_concordance_scatter <- function(df, x_col, y_col, 
                                     x_label, y_label,
                                     concordant_col) {
  
  df <- df %>%
    mutate(
      concordant = .data[[concordant_col]],
      direction  = case_when(
        .data[[x_col]] > 0 & .data[[y_col]] > 0 ~ "Higher in GABAergic (both)",
        .data[[x_col]] < 0 & .data[[y_col]] < 0 ~ "Higher in Glutamatergic (both)",
        TRUE                                      ~ "Discordant"
      ),
      direction = factor(direction, levels = c(
        "Higher in GABAergic (both)",
        "Higher in Glutamatergic (both)",
        "Discordant"
      ))
    )
  
  rho <- cor(df[[x_col]], df[[y_col]], method = "spearman", use = "complete.obs")
  
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_point(aes(color = direction), size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c(
        "Higher in GABAergic (both)"     = "#882255",
        "Higher in Glutamatergic (both)" = "#009988",
        "Discordant"                     = "grey70"
      ),
      name = NULL
    ) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.1, vjust = 1.5,
             label = paste0("Spearman ρ = ", round(rho, 2)),
             size = 2.8) +
    labs(
      x = paste0("Z-score diff (", x_label, ")"),
      y = paste0("Z-score diff (", y_label, ")")
    ) +
    theme_classic() +
    theme(
      axis.text       = element_text(size = 7),
      axis.title      = element_text(size = 8),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      aspect.ratio    = 1
    )
}

p_ss_vs_aging <- make_concordance_scatter(
  df            = df_concordance,
  x_col         = "diff_z_smartseq",
  y_col         = "diff_z_aging",
  x_label       = "Smart-seq",
  y_label       = "Jin et al. - Aging",
  concordant_col = "concordant_aging"
)

p_ss_vs_stress <- make_concordance_scatter(
  df            = df_concordance,
  x_col         = "diff_z_smartseq",
  y_col         = "diff_z_stress",
  x_label       = "Smart-seq",
  y_label       = "Hing et al. - Stress",
  concordant_col = "concordant_stress"
)

combined_concordance <- p_ss_vs_aging + p_ss_vs_stress +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

print(combined_concordance)

ggsave(
  "Figuras/Figure_S8AB_New.svg",
  plot   = combined_concordance,
  width  = 7,
  height = 3.5,
  units  = "in"
)

supp_table = df_concordance[,c(1,2,4,5,6,7)]
write.csv(supp_table,"data/Supplementary_Data_4_New.csv")


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

RP <- read.csv("data/RP.csv")
RP <- as.character(RP$x)

inhibitory <- c("Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip")

genes_of_interest <- c("Rps27", "Rpl21")

compute_specificity_scores <- function(seurat_obj, rp_genes, group_by = "predicted.id") {
  
  genes_present <- intersect(rp_genes, rownames(seurat_obj))
  cat(paste("Genes presentes:", length(genes_present), "de", length(rp_genes), "\n"))
  
  avg <- AverageExpression(
    object   = seurat_obj,
    features = genes_present,
    group.by = group_by,
    assay    = "RNA"
  )$RNA
  
  avg       <- avg[, colSums(avg) > 0, drop = FALSE]
  avg_dense <- as.matrix(avg)
  
  df        <- as.data.frame(avg_dense)
  df$mean_all <- rowMeans(df, na.rm = TRUE)
  
  subtypes     <- colnames(avg_dense)
  results_list <- list()
  
  for (subtype in subtypes) {
    model_data <- data.frame(
      Average_All  = df$mean_all,
      Mean_Subtype = df[[subtype]],
      row.names    = rownames(df)
    )
    model_data_clean <- na.omit(model_data)
    if (nrow(model_data_clean) > 1) {
      model    <- lm(Mean_Subtype ~ Average_All, data = model_data_clean)
      temp_res <- rep(NA_real_, nrow(df))
      names(temp_res) <- rownames(df)
      temp_res[rownames(model_data_clean)] <- residuals(model)
      results_list[[subtype]] <- temp_res
    }
  }
  
  residuals_df <- do.call(rbind, lapply(results_list, function(x) {
    matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }))
  rownames(residuals_df) <- names(results_list)
  
  global_mean <- mean(as.matrix(residuals_df), na.rm = TRUE)
  global_sd   <- sd(as.matrix(residuals_df),   na.rm = TRUE)
  
  specificity_scores <- (residuals_df - global_mean) / global_sd
  
  return(specificity_scores)
}

RP_smartseq <- gsub("Rack1", "Gnb2l1", RP)

load("data/Smartseq2dataset_seurat_filtered.RData")
DefaultAssay(seurat_merged) <- "RNA"
cat("Construyendo matriz Smart-seq (lógica Fig 3B)...\n")

subclass_pass <- list(
  c("L5 PT CTX","L5 IT CTX","L4/5 IT CTX","L6 IT CTX","L6 CT CTX","L5 NP CTX","Pvalb","Vip","L2/3 IT CTX","Lamp5","Sst","Sst Chodl","Sncg","Car3","L6b CTX"),
  c("L4/5 IT CTX","L2/3 IT CTX","Lamp5","L5 IT CTX","L6 IT CTX","Vip","L6 CT CTX","L5 PT CTX","Car3","L5 NP CTX"),
  c("Car3","L2/3 IT CTX"),
  c("Vip","L5 NP CTX","L4/5 IT CTX","Sst","L6 IT CTX","L5 IT CTX","Lamp5","L2/3 IT CTX","L6 CT CTX","L6b CTX"),
  c("L5 PT CTX","L4 RSP-ACA"),
  c("L2/3 IT CTX","L6 IT CTX"),
  c("Sst","L5 NP CTX","L5 IT CTX","Vip","Lamp5","L4/5 IT CTX","Pvalb","L2/3 IT CTX","L6 IT CTX","L6 CT CTX","L6b CTX","Sncg"),
  c("DG","CA1-ProS","CA3","Lamp5","Vip","Sncg","Pvalb","Sst"),
  c("Sst","Sst Chodl","Sncg","Vip","L6 CT CTX","L4/5 IT CTX","L5 IT CTX","L2/3 IT CTX","Lamp5","L5 NP CTX","Pvalb","L6 IT CTX","L5 PT CTX"),
  c("L2/3 IT CTX","L4/5 IT CTX"),
  c("L4/5 IT CTX","L2/3 IT CTX","L5 IT CTX","L2/3 IT PPP","L6 CT CTX","L5 NP CTX","L6 IT CTX","L6b CTX","Sst","Lamp5","Vip","Pvalb")
)
regions <- c("VISp","VIS","SSs","SSp","RSPv","RSP","MOp","HIP","ALM","AI","ACA")

seurat_split <- SplitObject(seurat_merged, split.by = "region_label")
seurat_split <- lapply(seurat_split, function(x) { DefaultAssay(x) <- "RNA"; x })

m <- matrix(nrow = length(RP_smartseq), ncol = 88)
colnames(m) <- unlist(lapply(seq_along(subclass_pass), function(i) paste(subclass_pass[[i]], regions[i], sep = "_")))
rownames(m) <- RP_smartseq

z <- 0
for (i in seq(1, 11)) {
  metadata_i <- as.data.frame(seurat_split[[i]]@meta.data)
  data_i     <- as.data.frame(GetAssayData(object = seurat_split[[i]], layer = "data"))
  for (p in seq_len(length(subclass_pass[[i]]))) {
    z <- z + 1
    cells      <- metadata_i[metadata_i$subclass_label %in% subclass_pass[[i]][p], ]
    datacellRP <- data_i[RP_smartseq, colnames(data_i) %in% rownames(cells)]
    m[, z]     <- rowMeans(datacellRP)
  }
}
remove(seurat_merged, seurat_split); gc()

rownames(m) <- gsub("Gnb2l1", "Rack1", rownames(m))

subtypes_ss <- unique(gsub("_.*", "", colnames(m)))
m_ss <- sapply(subtypes_ss, function(st) {
  cols <- grep(paste0("^", st, "_"), colnames(m), value = TRUE)
  rowMeans(m[, cols, drop = FALSE], na.rm = TRUE)
})

df_ss          <- as.data.frame(m_ss)
df_ss$mean_all <- rowMeans(as.data.frame(m), na.rm = TRUE)  # igual que Fig 3B: mean sobre las 88 cols

results_list_ss <- list()
for (subtype in subtypes_ss) {
  model_data <- data.frame(
    Average_All  = df_ss$mean_all,
    Mean_Subtype = df_ss[[subtype]],
    row.names    = rownames(df_ss)
  )
  model_data_clean <- na.omit(model_data)
  if (nrow(model_data_clean) > 1) {
    model    <- lm(Mean_Subtype ~ Average_All, data = model_data_clean)
    temp_res <- rep(NA_real_, nrow(df_ss))
    names(temp_res) <- rownames(df_ss)
    temp_res[rownames(model_data_clean)] <- residuals(model)
    results_list_ss[[subtype]] <- temp_res
  }
}

residuals_ss <- do.call(rbind, lapply(results_list_ss, function(x) {
  matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
}))
rownames(residuals_ss) <- names(results_list_ss)

global_mean_ss <- mean(as.matrix(residuals_ss), na.rm = TRUE)
global_sd_ss   <- sd(as.matrix(residuals_ss),   na.rm = TRUE)
scores_smartseq <- (residuals_ss - global_mean_ss) / global_sd_ss


# Jin Adult
seurat_adult <- readRDS("data/Normalized_Classified_Seurat_Adult")
meta_paper   <- data.table::fread("data/metadata_aging.csv")
metadata_adult <- seurat_adult@meta.data
metadata_adult$library_prep <- sapply(rownames(metadata_adult), function(x) strsplit(x, "-")[[1]][2])
metadata_adult$barcode <- rownames(metadata_adult)
merged_adult <- merge(metadata_adult, meta_paper, by = "library_prep", all.x = TRUE)
rownames(merged_adult) <- merged_adult$barcode
seurat_adult@meta.data <- merged_adult[rownames(seurat_adult@meta.data), ]
cat("Calculando scores Jin Adult...\n")
scores_adult <- compute_specificity_scores(seurat_adult, RP)
remove(seurat_adult); gc()

# Jin Aged
seurat_aged <- readRDS("data/Normalized_Classified_Seurat_Aged")
metadata_aged <- seurat_aged@meta.data
metadata_aged$library_prep <- sapply(rownames(metadata_aged), function(x) strsplit(x, "-")[[1]][2])
metadata_aged$barcode <- rownames(metadata_aged)
merged_aged <- merge(metadata_aged, meta_paper, by = "library_prep", all.x = TRUE)
rownames(merged_aged) <- merged_aged$barcode
seurat_aged@meta.data <- merged_aged[rownames(seurat_aged@meta.data), ]
cat("Calculando scores Jin Aged...\n")
scores_aged <- compute_specificity_scores(seurat_aged, RP)
remove(seurat_aged); gc()

# Hing Stress
seurat_stress <- LoadSeuratRds("data/Hing_2024_PFC_stress/SeuratProject.Rds")
cat("Calculando scores Hing Stress...\n")
scores_stress <- compute_specificity_scores(seurat_stress, RP)
remove(seurat_stress); gc()

extract_scores_long <- function(scores_mat, dataset_label) {
  as.data.frame(scores_mat) %>%
    rownames_to_column("Subtype") %>%
    pivot_longer(-Subtype, names_to = "Gene", values_to = "Score") %>%
    mutate(Dataset = dataset_label)
}

df_smartseq <- extract_scores_long(scores_smartseq, "Smart-seq")

df_all <- bind_rows(
  df_smartseq,
  extract_scores_long(scores_adult,  "Jin et al. - Adult"),
  extract_scores_long(scores_aged,   "Jin et al. - Aged"),
  extract_scores_long(scores_stress, "Hing et al. - Stress")
) %>%
  filter(Gene %in% genes_of_interest) %>%
  mutate(
    Class         = ifelse(Subtype %in% inhibitory, "GABAergic", "Glutamatergic"),
    Dataset       = factor(Dataset, levels = c("Smart-seq",
                                               "Jin et al. - Adult",
                                               "Jin et al. - Aged",
                                               "Hing et al. - Stress")),
    Score_clamped = pmax(pmin(Score, 4), -4)
  )

subtype_order <- df_all %>%
  filter(Gene == "Rps27") %>%
  group_by(Subtype, Class) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  arrange(Class, desc(mean_score)) %>%
  pull(Subtype) %>%
  unique()

n_inhib <- sum(subtype_order %in% inhibitory)

df_all <- df_all %>%
  mutate(Subtype = factor(Subtype, levels = subtype_order))

p_heatmap <- ggplot(df_all, aes(x = Dataset, y = Subtype, fill = Score_clamped)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(
    data = df_all %>% filter(abs(Score) > 2.5),
    aes(label = round(Score, 1)),
    size = 2.5, color = "white", fontface = "bold"
  ) +
  geom_hline(yintercept = n_inhib + 0.5,
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("#009988", "white", "#882255"),
    limits = c(-4, 4),
    name   = "Specificity\nScore"
  ) +
  facet_wrap(~ Gene, ncol = 2) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 8),
    strip.text      = element_text(size = 10, face = "bold"),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    panel.grid      = element_blank()
  )

print(p_heatmap)

ggsave(
  "Figuras/Figure_S8C_NEW.svg",
  plot   = p_heatmap,
  width  = 8,
  height = 6,
  units  = "in"
)

