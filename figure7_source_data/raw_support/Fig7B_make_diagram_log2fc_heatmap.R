#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(RColorBrewer)
  library(scales)
})

base_dir <- "/Users/xq/Desktop/睾丸文章总/最新内容汇总/figure汇总/figure7"
in_path <- file.path(base_dir, "late阶段_代谢通路图基因表达_提取.csv")
out_dir <- file.path(base_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(in_path, show_col_types = FALSE)

gene_levels <- df %>%
  arrange(diagram_order) %>%
  pull(gene)

df <- df %>%
  mutate(
    pathway = factor(pathway, levels = c("Glycolysis", "TCA cycle", "Oxidative phosphorylation")),
    comparison = "Disease vs Ctrl",
    gene = factor(gene, levels = rev(gene_levels)),
    log2fc_clip = squish(log2fc, range = c(-2, 2))
  )

palette <- rev(brewer.pal(11, "Spectral"))

p <- ggplot(df, aes(x = comparison, y = gene, fill = log2fc_clip)) +
  geom_tile(width = 0.95, height = 0.95, color = "white", linewidth = 0.15) +
  facet_grid(pathway ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = palette,
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    name = expression(paste("Gene expression\n", Log[2], "FC"))
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      direction = "horizontal",
      barwidth = unit(3.5, "in"),
      barheight = unit(0.18, "in"),
      ticks = TRUE
    )
  ) +
  labs(x = NULL, y = NULL, title = "Late primary SPCs: diagram genes (Disease vs Ctrl)") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0)
  )

ggsave(
  filename = file.path(out_dir, "late_stage_diagram_genes_log2fc_heatmap.pdf"),
  plot = p,
  width = 7,
  height = 14,
  useDingbats = FALSE
)
ggsave(
  filename = file.path(out_dir, "late_stage_diagram_genes_log2fc_heatmap.png"),
  plot = p,
  width = 7,
  height = 14,
  dpi = 300
)

message("Done. Outputs in: ", out_dir)

