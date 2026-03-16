#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path))
} else {
  normalizePath(getwd())
}

data_dir <- normalizePath(file.path(script_dir, "..", "fig1_source_data"))
out_dir <- file.path(script_dir, "fig1_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(plot_obj, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot_obj, width = width, height = height, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot_obj, width = width, height = height, dpi = 300)
}

group_pref <- c("Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS")
celltype_pref <- c("SSCs", "SPGs", "SPCs", "STs", "Spermatids", "LCs", "ECs", "Myeloid", "Lym", "Myoid")
marker_pref <- c(
  "VIM", "SOX9", "VWF", "DLK1", "PTPRC", "NOTCH3",
  "DDX4", "FGFR3", "KIT", "SYCP3", "TEX29", "PRM3",
  "CPEB1", "MAGEB2", "DMRT1", "SOHLH1", "SOHLH2", "BRCA1", "DMC1",
  "MLH3", "PSMC3IP", "WT1", "CLDN11", "AMH", "S100A13", "ENO1", "BEX1",
  "CST9L", "HOPX", "SPAG6", "ACTL7B", "ATG4D", "ZPBP", "DNAH6", "SPATA18",
  "HOOK1", "TNP1", "NR2F2", "CFD", "CLDN5", "CDH5", "PECAM1", "LYZ",
  "CD14", "MRC1", "C1QC", "AIF1", "CD3D", "CD79A", "ACTA2", "MYH11"
)

# Panel B
umap_group <- read_csv(file.path(data_dir, "Fig1B_UMAP_by_group.csv"), show_col_types = FALSE) %>%
  filter(!is.na(UMAP1), !is.na(UMAP2), !is.na(group))

group_levels <- c(group_pref[group_pref %in% unique(umap_group$group)], setdiff(unique(umap_group$group), group_pref))
umap_group <- umap_group %>% mutate(group = factor(group, levels = group_levels))

p_b <- ggplot(umap_group, aes(UMAP1, UMAP2, color = group)) +
  geom_point(size = 0.08, alpha = 0.65) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "right")
save_panel(p_b, "Fig1B_UMAP_by_group", 5, 4)

# Panel C
umap_celltype <- read_csv(file.path(data_dir, "Fig1C_UMAP_by_celltype.csv"), show_col_types = FALSE) %>%
  filter(!is.na(UMAP1), !is.na(UMAP2), !is.na(cell_type))

ct_levels <- c(celltype_pref[celltype_pref %in% unique(umap_celltype$cell_type)], setdiff(unique(umap_celltype$cell_type), celltype_pref))
umap_celltype <- umap_celltype %>% mutate(cell_type = factor(cell_type, levels = ct_levels))

p_c <- ggplot(umap_celltype, aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.08, alpha = 0.65) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "right")
save_panel(p_c, "Fig1C_UMAP_by_celltype", 5.5, 4.5)

# Panel D
set.seed(1)
feature <- read_csv(file.path(data_dir, "Fig1D_marker_feature_long.csv"), show_col_types = FALSE) %>%
  filter(!is.na(UMAP1), !is.na(UMAP2), !is.na(gene), !is.na(expression)) %>%
  group_by(gene) %>%
  slice_sample(n = min(n(), 40000)) %>%
  ungroup()

marker_levels <- c(marker_pref[marker_pref %in% unique(feature$gene)], setdiff(unique(feature$gene), marker_pref))
feature <- feature %>% mutate(gene = factor(gene, levels = marker_levels))

p_d <- ggplot(feature, aes(UMAP1, UMAP2, color = expression)) +
  geom_point(size = 0.05, alpha = 0.75) +
  facet_wrap(~gene, ncol = 4) +
  coord_equal() +
  scale_color_gradientn(colors = c("#d9d9d9", "#fdae6b", "#d7301f"), trans = "sqrt") +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 9)
  )
save_panel(p_d, "Fig1D_featureplots", 10, 7)

# Panel E
comp <- read_csv(file.path(data_dir, "Fig1E_group_composition_percent_long.csv"), show_col_types = FALSE) %>%
  filter(!is.na(group), !is.na(celltype), !is.na(percent)) %>%
  mutate(
    percent_num = parse_number(as.character(percent)),
    percent_num = ifelse(grepl("%", as.character(percent)), percent_num / 100, percent_num)
  ) %>%
  filter(group != "All samples")

comp_levels <- c(group_pref[group_pref %in% unique(comp$group)], setdiff(unique(comp$group), group_pref))
ct2_levels <- c(celltype_pref[celltype_pref %in% unique(comp$celltype)], setdiff(unique(comp$celltype), celltype_pref))
comp <- comp %>%
  mutate(
    group = factor(group, levels = comp_levels),
    celltype = factor(celltype, levels = ct2_levels)
  )

p_e <- ggplot(comp, aes(x = group, y = percent_num, fill = celltype)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.15) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Fraction", fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(p_e, "Fig1E_composition", 5, 4)

# Panel G
dot_path_full <- file.path(data_dir, "Fig1G_dotplot_full_genes.csv")
dot_path_fallback <- file.path(data_dir, "Fig1G_dotplot_available_genes.csv")
dot <- if (file.exists(dot_path_full)) {
  read_csv(dot_path_full, show_col_types = FALSE)
} else {
  read_csv(dot_path_fallback, show_col_types = FALSE)
}

dot <- dot %>%
  filter(!is.na(cell_type), !is.na(gene), !is.na(avg_expr), !is.na(pct_expr))
if (max(dot$pct_expr, na.rm = TRUE) > 1) {
  dot <- dot %>% mutate(pct_expr = pct_expr / 100)
}

ct_dot_levels <- c(celltype_pref[celltype_pref %in% unique(dot$cell_type)], setdiff(unique(dot$cell_type), celltype_pref))
gene_levels <- c(marker_pref[marker_pref %in% unique(dot$gene)], setdiff(unique(dot$gene), marker_pref))
dot <- dot %>%
  mutate(
    cell_type = factor(cell_type, levels = ct_dot_levels),
    gene = factor(gene, levels = gene_levels)
  )

p_g <- ggplot(dot, aes(x = gene, y = cell_type)) +
  geom_point(aes(size = pct_expr, color = avg_expr), alpha = 0.95) +
  scale_size(range = c(0.4, 6.2), limits = c(0, 1), labels = percent_format(accuracy = 1), name = "% expressed") +
  scale_color_gradientn(colors = c("#2166ac", "#f7f7f7", "#b2182b"), name = "Avg expression") +
  theme_bw(base_size = 10) +
  labs(x = "Gene", y = "Cell type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank()
  )
save_panel(p_g, "Fig1G_dotplot", 11, 5)

final <- (p_b | p_c) / p_d / (p_e | p_g) +
  plot_layout(heights = c(1, 1.35, 1.1), widths = c(1, 1.4))
save_panel(final, "Figure1_rebuild_from_source_data", 14, 15)

message("Fig1 plots saved to ", out_dir)
