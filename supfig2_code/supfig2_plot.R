#!/usr/bin/env Rscript

## Supplementary Figure 2 plotting script (Panels A–B)
## Uses source data in ../supfig2_source_data
## Output: supfig2_plots/ (panel PNG/PDF + combined figure)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(cowplot)
  library(patchwork)
  library(scales)
})

# ---- paths ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg) else NULL
script_dir <- if (!is.null(script_path)) dirname(normalizePath(script_path)) else getwd()

data_dir <- normalizePath(file.path(script_dir, "..", "supfig2_source_data"))
out_dir <- file.path(script_dir, "supfig2_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Panel A: feature plots ----
feature_path <- file.path(data_dir, "SupFig2A_featureplot_cells.csv")
feature <- read_csv(feature_path, show_col_types = FALSE)

gene_order_A <- c("SYCP3", "PIWIL1", "TNP1", "PRM2")
group_order_A <- c("Ctrl", "OA", "KS", "iNOA_B", "AZFc_Del")

feature <- feature %>%
  filter(group %in% group_order_A, gene %in% gene_order_A) %>%
  mutate(
    gene = factor(gene, levels = gene_order_A),
    group = factor(group, levels = group_order_A)
  )

pA <- ggplot(feature, aes(x = UMAP1, y = UMAP2, color = expr)) +
  geom_point(size = 0.18, alpha = 0.9) +
  scale_color_gradient(
    low = "#f2f2f2",
    high = "#c61b1b",
    limits = c(0, 5),
    oob = scales::squish,
    na.value = "#f2f2f2"
  ) +
  facet_grid(rows = vars(group), cols = vars(gene)) +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 9),
    panel.spacing = unit(0.25, "lines")
  )

# ---- Panel B: dotplot ----
dot_path <- file.path(data_dir, "SupFig2B_dotplot_groups.csv")
dot <- read_csv(dot_path, show_col_types = FALSE)

gene_order_B <- c("PRM2", "TNP1", "PIWIL1", "SYCP3")
group_order_B <- c("Ctrl", "iNOA_B", "OA", "AZFc_Del", "KS")

dot <- dot %>%
  mutate(
    gene = factor(gene, levels = gene_order_B),
    group = factor(group, levels = group_order_B)
  )

pB <- ggplot(dot, aes(x = group, y = gene)) +
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", limits = c(0, 1)) +
  scale_size(range = c(1, 6), limits = c(0, 100), breaks = c(20, 40, 60, 80, 100)) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(size = 0.2)
  )

# ---- Combine ----
combined <- plot_grid(
  pA, pB,
  labels = c("A", "B"),
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(3.2, 1.1)
)

# ---- Save outputs ----
ggsave(file.path(out_dir, "SupFig2A.png"), pA, width = 8.5, height = 7.5, dpi = 300)
ggsave(file.path(out_dir, "SupFig2B.png"), pB, width = 4.2, height = 2.2, dpi = 300)
ggsave(file.path(out_dir, "SupFig2_combined.png"), combined, width = 8.5, height = 9.2, dpi = 300)

ggsave(file.path(out_dir, "SupFig2A.pdf"), pA, width = 8.5, height = 7.5)
ggsave(file.path(out_dir, "SupFig2B.pdf"), pB, width = 4.2, height = 2.2)
ggsave(file.path(out_dir, "SupFig2_combined.pdf"), combined, width = 8.5, height = 9.2)
