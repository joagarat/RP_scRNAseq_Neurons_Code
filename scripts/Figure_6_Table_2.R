
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggbeeswarm)
library(tibble)
library(scales)
library(gt)
library(janitor) 

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

load_and_process_data <- function(base_path, rp_genes) {
  
  cat(paste("Procesando datos desde:", base_path, "\n"))
  
  filelist <- list.files(path = base_path, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  
  if (length(filelist) == 0) {
    stop("No se encontraron archivos .csv en la ruta especificada: ", base_path)
  }
  
  regions <- sapply(stringr::str_split(filelist, "/"), function(path) path[length(path) - 1])
  
  df_input_list <- lapply(seq_along(filelist), function(i) {
    df <- data.table::fread(filelist[i])
    df$region <- regions[i]
    df$comparison_raw <- gsub("^.*/|\\.csv$", "", filelist[i])
    df
  })
  
  df_allregions <- dplyr::bind_rows(df_input_list)
  
  df_allregions_rp <- df_allregions[df_allregions$V1 %in% rp_genes, ]
  df_allregions_rp$log2FoldChange <- df_allregions_rp$log2FoldChange * -1 
  
  df_corrected <- df_allregions_rp %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      standardized = list(standardize_comparison(comparison_raw, log2FoldChange)),
      comparison_id = standardized$name,
      log2FoldChange = standardized$log2fc 
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(V1, comparison_id, region, log2FoldChange, padj) 
  
  return(df_corrected)
}

RP_df <- read.csv("data/RP.csv")
RP <- as.character(RP_df$x)
RP <- gsub("Rack1", "Gnb2l1", x = RP)
path_10x <- "data/Chromium_DEGs/"
path_smartseq <- "data/Smartseq_DEGs/"
df_corrected_10x <- load_and_process_data(base_path = path_10x, rp_genes = RP)
df_corrected_smartseq <- load_and_process_data(base_path = path_smartseq, rp_genes = RP)
df_corrected_10x <- df_corrected_10x %>%
  mutate(V1 = ifelse(V1 == "Rack1", "Gnb2l1", V1))
inhibitory_subclasses <- c("Sst", "Pvalb", "Vip", "Lamp5", "Sncg", "Sst Chodl")
excitatory_subclasses <- c("CA1-ProS","CA3", "Car3", "DG","L2_3 IT CTX", "L2_3 IT PPP","L4 RSP-ACA",
                           "L4_5 IT CTX", "L5 IT CTX", "L5 NP CTX","L5 PT CTX","L6 CT CTX", 
                           "L6 IT CTX" , "L6b CTX")
all_subclasses_list <- unique(c(inhibitory_subclasses, excitatory_subclasses)) # Usamos unique por si acaso
subclass_info <- data.frame(subclass = all_subclasses_list) %>%
  mutate(
    Class = ifelse(subclass %in% inhibitory_subclasses, "Inhibitory", "Excitatory")
  )
rownames(subclass_info) <- subclass_info$subclass
concordance_data <- inner_join(
  df_corrected_smartseq,
  df_corrected_10x,
  by = c("V1", "comparison_id", "region"),
  suffix = c("_smartseq", "_10x")
) %>%
  drop_na(log2FoldChange_smartseq, padj_smartseq, log2FoldChange_10x, padj_10x) %>%
  rename(gene = V1)
logfc_threshold <- 0.58496250072
padj_threshold <- 0.05

scatter_data <- concordance_data %>%
  mutate(
    concordance_type = case_when(
      padj_smartseq < padj_threshold & abs(log2FoldChange_smartseq) > logfc_threshold &
        padj_10x < padj_threshold & abs(log2FoldChange_10x) > logfc_threshold &
        sign(log2FoldChange_smartseq) == sign(log2FoldChange_10x) ~ "Concordant RP",
      padj_smartseq < padj_threshold & abs(log2FoldChange_smartseq) > logfc_threshold &
        padj_10x < padj_threshold & abs(log2FoldChange_10x) > logfc_threshold &
        sign(log2FoldChange_smartseq) != sign(log2FoldChange_10x) ~ "Discordant RP",
      TRUE ~ "Non-Replicated RP"
    ),
    concordance_type = factor(concordance_type, levels = c("Non-Replicated RP", "Discordant RP", "Concordant RP"))
  ) %>%
  arrange(concordance_type)
correlation <- cor(
  scatter_data$log2FoldChange_smartseq,
  scatter_data$log2FoldChange_10x,
  method = "spearman"
)

panel_a_scatterplot <- ggplot(scatter_data, aes(x = log2FoldChange_smartseq, y = log2FoldChange_10x)) +
  geom_point(aes(color = concordance_type), alpha = 0.7, size = 2) +
  scale_color_manual(
    name = "RP Concordance", 
    values = c(
      "Concordant RP" = "#e41a1c",       
      "Discordant RP" = "#377eb8",      
      "Non-Replicated RP" = "grey80"    
    )
  ) +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  
  annotate(
    "text", 
    x = min(scatter_data$log2FoldChange_smartseq, na.rm = TRUE), 
    y = max(scatter_data$log2FoldChange_10x, na.rm = TRUE), 
    label = paste("Spearman's ρ =", round(correlation, 2)), 
    hjust = 0, vjust = 1, 
    size = 3 
  ) +
  
  theme_classic() + 
  theme(
    legend.position = "right", 
    aspect.ratio = 1,
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) + 
  
  labs(
    title = NULL,
    x = "Subclass-Level log2FC (Smart-seq)",
    y = "Subclass-Level log2FC (10x Genomics)"
  )

print(panel_a_scatterplot)

ggsave(
  "Figuras/Figura_6/Figure_6A.svg", 
  plot = panel_a_scatterplot, 
  width = 4.5,  
  height = 3.5, 
  units = "in"
)

concordant_degs <- scatter_data %>%
  filter(concordance_type == "Concordant RP")
winner_loser_data <- concordant_degs %>%
  separate(comparison_id, into = c("subclass1", "subclass2"), sep = "vs", remove = FALSE) %>%
  mutate(
    winner_subclass = case_when(
      log2FoldChange_smartseq < 0 ~ subclass2,
      log2FoldChange_smartseq > 0 ~ subclass1
    ),
    loser_subclass = case_when(
      log2FoldChange_smartseq < 0 ~ subclass1,
      log2FoldChange_smartseq > 0 ~ subclass2
    )
  )
dot_plot_winners <- winner_loser_data %>%
  group_by(gene, winner_subclass) %>%
  summarise(
    n_concordant = n(),
    avg_abs_log2FC = mean(abs(log2FoldChange_smartseq), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(subclass = winner_subclass)

dot_plot_losers <- winner_loser_data %>%
  group_by(gene, loser_subclass) %>%
  summarise(
    n_concordant = n(),
    avg_abs_log2FC = mean(abs(log2FoldChange_smartseq), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(subclass = loser_subclass)

gene_activity <- bind_rows(
  dot_plot_winners %>% select(gene, n_concordant),
  dot_plot_losers %>% select(gene, n_concordant)
) %>%
  group_by(gene) %>%
  summarise(total_activity = sum(n_concordant)) %>%
  arrange(desc(total_activity))

gene_order <- gene_activity %>%
  slice_head(n = 20) %>%
  pull(gene)

subclass_activity_tibble <- bind_rows(
  dot_plot_winners %>% rename(count = n_concordant),
  dot_plot_losers %>% rename(count = n_concordant)
) %>%
  group_by(subclass, gene) %>%
  summarise(total_activity = sum(count), .groups = "drop") %>%
  pivot_wider(id_cols = subclass, names_from = gene, values_from = total_activity, values_fill = 0)

subclass_activity_matrix <- subclass_activity_tibble %>%
  column_to_rownames(var = "subclass") %>%
  as.matrix()
subclass_order <- hclust(dist(subclass_activity_matrix))$order
subclass_levels <- rownames(subclass_activity_matrix)[subclass_order]
dot_plot_winners$gene <- factor(dot_plot_winners$gene, levels = rev(gene_order))
dot_plot_winners$subclass <- factor(dot_plot_winners$subclass, levels = subclass_levels)
dot_plot_losers$gene <- factor(dot_plot_losers$gene, levels = rev(gene_order))
dot_plot_losers$subclass <- factor(dot_plot_losers$subclass, levels = subclass_levels)


panel_b1_dotplot_winners <- ggplot(dot_plot_winners %>% filter(gene %in% gene_order), 
                                   aes(x = subclass, y = gene)) +
  geom_point(aes(size = n_concordant, color = avg_abs_log2FC)) +
  
  scale_color_gradientn(
    name = "Mean |log2FC|",
    colors = c("#FCAE91", "#99000D"),
    values = scales::rescale(c(0, 1.5)),
    limits = c(0, 1.5), oob = scales::squish
  ) +
  scale_size_continuous(
    name = "# Concordant OE", 
    limits = c(0, 50), range = c(1, 5)
  ) +
  
  labs(title = "Concordant Overexpression", x = NULL, y = NULL) +
  
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
  ) +
  
  guides(
    color = guide_colorbar(order = 1, barheight = unit(4, "lines"), barwidth = unit(0.5, "lines")),
    size = guide_legend(order = 2)
  )

panel_b2_dotplot_losers <- ggplot(dot_plot_losers %>% filter(gene %in% gene_order), 
                                  aes(x = subclass, y = gene)) +
  geom_point(aes(size = n_concordant, color = avg_abs_log2FC)) +
  
  scale_color_gradientn(
    name = "Mean |log2FC|",
    colors = c("#6BAED6", "#08306B"),
    values = scales::rescale(c(0, 1.5)),
    limits = c(0, 1.5), oob = scales::squish
  ) +
  scale_size_continuous(
    name = "# Concordant UE", 
    limits = c(0, 50), range = c(1, 5)
  ) +
  
  labs(title = "Concordant Underexpression", x = NULL, y = NULL) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_blank(),     # Sin eje Y (compartido)
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
  ) +
  
  guides(
    color = guide_colorbar(order = 1, barheight = unit(4, "lines"), barwidth = unit(0.5, "lines")),
    size = guide_legend(order = 2)
  )

panel_b_final <- panel_b1_dotplot_winners + panel_b2_dotplot_losers

panel_b_final <- panel_b_final + 
  plot_annotation(
    caption = "Neuronal Subclass",
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = 10, face = "bold") 
    )
  )

print(panel_b_final)

ggsave(
  "Figuras/Figura_6/Figure_6B_C.svg",
  plot = panel_b_final,
  width = 7.5,  
  height = 5.5,
  units = "in"
)

smartseq_subclass_keys <- df_corrected_smartseq %>%
  filter(comparison_id != "GABAergicvsGlutamatergic") %>%
  distinct(region, comparison_id)

x10x_subclass_keys <- df_corrected_10x %>%
  filter(comparison_id != "GABAergicvsGlutamatergic") %>%
  distinct(region, comparison_id)

shared_subclass_comparisons <- inner_join(
  smartseq_subclass_keys,
  x10x_subclass_keys,
  by = c("region", "comparison_id")
)

padj_threshold <- 0.05
logfc_threshold <- 0.58496250072

smartseq_shared_findings <- df_corrected_smartseq %>%
  semi_join(shared_subclass_comparisons, by = c("region", "comparison_id")) %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > logfc_threshold) %>%
  group_by(region) %>%
  summarise(`Subclass RP DEGs Smartseq` = n(), .groups = 'drop')

x10x_shared_findings <- df_corrected_10x %>%
  semi_join(shared_subclass_comparisons, by = c("region", "comparison_id")) %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > logfc_threshold) %>%
  group_by(region) %>%
  summarise(`Subclass RP DEGs 10x Genomics` = n(), .groups = 'drop')

concordance_counts <- scatter_data %>%
  filter(comparison_id != "GABAergicvsGlutamatergic") %>%
  group_by(region, concordance_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = concordance_type, values_from = count, values_fill = 0) %>%
  select(
    region,
    `Shared Congruent` = `Concordant RP`,
    `Incongruent` = `Discordant RP`
  )

final_concordance_table_data <- concordance_counts %>%
  full_join(smartseq_shared_findings, by = "region") %>%
  full_join(x10x_shared_findings, by = "region") %>%
  filter(region %in% shared_subclass_comparisons$region) %>%
  select(
    region,
    `Subclass RP DEGs Smartseq`,
    `Subclass RP DEGs 10x Genomics`,
    `Shared Congruent`,
    `Incongruent`
  ) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  arrange(region)

final_concordance_table_data_with_total <- final_concordance_table_data %>%
  adorn_totals("row")

print(final_concordance_table_data_with_total)


publication_concordance_table <- final_concordance_table_data_with_total %>%
  gt(rowname_col = "region") %>%
  tab_header(
    title = "Concordance of Subclass-Level DE RPs between Smart-seq and 10x Genomics Datasets"
  ) %>%
  cols_label(
    `Subclass RP DEGs Smartseq` = "Total DE RPs (Smart-seq)",
    `Subclass RP DEGs 10x Genomics` = "Total DE RPs (10x)",
    `Shared Congruent` = "Shared & Concordant",
    `Incongruent` = "Shared & Discordant"
  ) %>%
  sub_missing(missing_text = "") %>%
  cols_align(
    align = "right",
    columns = where(is.numeric)
  ) %>%
  tab_style(
    style = list(cell_fill(color = "grey95"), cell_text(weight = "bold")),
    locations = list(
      cells_stub(rows = "Total"),
      cells_body(rows = "Total")
    )
  ) %>%
  tab_source_note(
    md("DE RPs were defined by *p*~adj~ < 0.05 and an absolute fold-change > 1.5.")
  ) %>%
  tab_options(
    table.border.top.width = px(3),
    table.border.bottom.width = px(3),
    column_labels.border.bottom.width = px(2)
  )

publication_concordance_table



gtsave(publication_concordance_table, "Figuras/Table2_Concordance_DEGs.png", zoom = 4)







