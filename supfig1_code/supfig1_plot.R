#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path))
} else {
  normalizePath(getwd())
}

data_dir <- normalizePath(file.path(script_dir, "..", "supfig1_source_data"))
out_dir <- file.path(script_dir, "supfig1_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(plot_obj, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot_obj, width = width, height = height, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot_obj, width = width, height = height, dpi = 300)
}

# Sup Fig 1A
set.seed(1)
feat <- read_csv(file.path(data_dir, "SupFig1A_featureplot_long.csv"), show_col_types = FALSE) %>%
  filter(!is.na(UMAP1), !is.na(UMAP2), !is.na(gene), !is.na(expression)) %>%
  group_by(gene) %>%
  slice_sample(n = min(n(), 35000)) %>%
  ungroup()

gene_levels <- unique(feat$gene)
feat <- feat %>% mutate(gene = factor(gene, levels = gene_levels))

p_a <- ggplot(feat, aes(UMAP1, UMAP2, color = expression)) +
  geom_point(size = 0.05, alpha = 0.75) +
  facet_wrap(~gene, ncol = 5) +
  coord_equal() +
  scale_color_gradientn(colors = c("#d9d9d9", "#fdae6b", "#d7301f"), trans = "sqrt") +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "right"
  )
save_panel(p_a, "SupFig1A_featureplots", 14, 10)

# Sup Fig 1B
heat <- read_csv(file.path(data_dir, "SupFig1B_heatmap_scaled_long.csv"), show_col_types = FALSE) %>%
  filter(!is.na(gene), !is.na(cell_type), !is.na(value))

celltype_pref <- c("SSCs", "SPGs", "SPCs", "STs", "Spermatids", "LCs", "ECs", "Myeloid", "Lym", "Myoid")
ct_levels <- c(celltype_pref[celltype_pref %in% unique(heat$cell_type)], setdiff(unique(heat$cell_type), celltype_pref))
gene_levels_hm <- unique(heat$gene)

heat <- heat %>%
  mutate(
    cell_type = factor(cell_type, levels = ct_levels),
    gene = factor(gene, levels = rev(gene_levels_hm))
  )

p_b <- ggplot(heat, aes(cell_type, gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Scaled\nexpr") +
  theme_minimal(base_size = 10) +
  labs(x = "Cell type", y = "Gene") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5),
    legend.position = "right"
  )
save_panel(p_b, "SupFig1B_heatmap_scaled", 9, 14)

combined <- p_a / p_b + plot_layout(heights = c(1.2, 1.05))
save_panel(combined, "SupFig1_combined", 14, 20)

message("SupFig1 plots saved to ", out_dir)
