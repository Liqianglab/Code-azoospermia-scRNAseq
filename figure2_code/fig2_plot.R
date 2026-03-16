#!/usr/bin/env Rscript

## Figure 2 plotting script (panels A–H)
## Public version: uses only source data in ../figure2_source_data.
## Output: fig2_plots/ (panel PDFs + combined Figure2 PDF/PNG)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(cowplot)
  library(patchwork)
  library(scales)
})

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && a != "") a else b

# ---- paths ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg) else NULL
script_dir <- if (!is.null(script_path)) dirname(normalizePath(script_path)) else getwd()

data_dir <- normalizePath(file.path(script_dir, "..", "figure2_source_data"))
out_dir <- file.path(script_dir, "fig2_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- colors ----
stage_order <- c("SSCs","SPGs","Early_primary_SPCs","Late_primary_SPCs",
                 "Round_Spermatids","Elongated_Spermatids","Sperm")
stage_colors <- c(
  SSCs = "#1b9e77",
  SPGs = "#d95f02",
  Early_primary_SPCs = "#7570b3",
  Late_primary_SPCs = "#e7298a",
  Round_Spermatids = "#66a61e",
  Elongated_Spermatids = "#e6ab02",
  Sperm = "#a6761d"
)

group_order <- c("Ctrl","iNOA_B","OA","AZFc_Del","KS")
group_colors <- c(
  Ctrl = "#1f77b4",
  iNOA_B = "#ff7f0e",
  OA = "#2ca02c",
  AZFc_Del = "#d62728",
  KS = "#9467bd"
)

tag_panel <- function(p, label) {
  ggdraw(p) + draw_label(label, x = 0, y = 1, hjust = -0.1, vjust = 1.2, fontface = "bold", size = 12)
}

# ---- Panel A: UMAP by stage + composition ----
umap_stage <- read_csv(file.path(data_dir, "Fig2A_UMAP_by_stage.csv"), show_col_types = FALSE)
umap_stage$stage <- factor(umap_stage$stage, levels = stage_order)
umap_stage$group <- factor(umap_stage$group, levels = group_order)

pA1 <- ggplot(umap_stage, aes(UMAP1, UMAP2, color = stage)) +
  geom_point(size = 0.2, alpha = 0.8) +
  scale_color_manual(values = stage_colors, drop = FALSE) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_blank())

comp_path <- file.path(data_dir, "Fig2A_group_stage_composition_percent_from_samples.csv")
comp <- read_csv(comp_path, show_col_types = FALSE)
comp$group <- factor(comp$group, levels = group_order)
comp$stage <- factor(comp$stage, levels = stage_order)

pA2 <- ggplot(comp, aes(x = group, y = percent, fill = stage)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = stage_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

pA <- pA1 + pA2 + plot_layout(widths = c(2.2, 1))

# ---- Panel B: UMAP by group ----
umap_group <- read_csv(file.path(data_dir, "Fig2B_UMAP_by_group.csv"), show_col_types = FALSE)
umap_group$group <- factor(umap_group$group, levels = group_order)

pB <- ggplot(umap_group, aes(UMAP1, UMAP2, color = group)) +
  geom_point(size = 0.2, alpha = 0.8) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_blank())

# ---- Panel C: dotplot ----
dot <- read_csv(file.path(data_dir, "Fig2C_dotplot_germline_markers.csv"), show_col_types = FALSE)
marker_file <- file.path(data_dir, "Fig2C_marker_order.csv")
if (file.exists(marker_file)) {
  markers <- read_csv(marker_file, show_col_types = FALSE) %>% pull(1) %>% as.character()
  markers <- markers[markers %in% unique(dot$gene)]
  if (length(markers) > 0) {
    dot$gene <- factor(dot$gene, levels = markers)
  }
}
dot$stage <- factor(dot$stage, levels = stage_order)

pC <- ggplot(dot, aes(x = gene, y = stage)) +
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_color_gradient(low = "#fde0dd", high = "#c51b8a") +
  scale_size(range = c(0.5, 4)) +
  theme_minimal(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(size = 0.2)
  )

# ---- Panel D: trajectory by group ----
traj_cells <- read_csv(file.path(data_dir, "Fig2D_pseudotime_cells.csv"), show_col_types = FALSE)
traj_cells$group <- factor(traj_cells$group, levels = group_order)
traj_cells$stage <- factor(traj_cells$stage, levels = stage_order)

traj_lines_path <- file.path(data_dir, "Fig2D_trajectory_lines.csv")
traj_lines <- if (file.exists(traj_lines_path)) read_csv(traj_lines_path, show_col_types = FALSE) else NULL

make_traj_panel <- function(grp) {
  df <- traj_cells %>% filter(group == grp)
  # find arrow start/end based on pseudotime quantiles
  q1 <- quantile(df$pseudotime, 0.1, na.rm = TRUE)
  q9 <- quantile(df$pseudotime, 0.9, na.rm = TRUE)
  start <- df %>% mutate(d = abs(pseudotime - q1)) %>% arrange(d) %>% slice(1)
  end <- df %>% mutate(d = abs(pseudotime - q9)) %>% arrange(d) %>% slice(1)

  p <- ggplot() +
    {if (!is.null(traj_lines)) geom_path(data = traj_lines, aes(x = X1, y = X2, group = ID), color = "grey70", size = 0.3)} +
    geom_point(data = df, aes(x = X1, y = X2, color = stage), size = 0.4, alpha = 0.9) +
    geom_segment(
      data = data.frame(x = start$X1, y = start$X2, xend = end$X1, yend = end$X2),
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.08, "inches")),
      color = "black"
    ) +
    scale_color_manual(values = stage_colors, drop = FALSE) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(grp)
  p
}

pD_list <- lapply(group_order, make_traj_panel)
pD <- wrap_plots(pD_list, nrow = 1)

# ---- Panel E: divergence along pseudotime ----
curve <- read_csv(file.path(data_dir, "Fig2E_entropy_curve.csv"), show_col_types = FALSE)
stage_bounds <- read_csv(file.path(data_dir, "Fig2E_stage_pseudotime_q10_q90.csv"), show_col_types = FALSE)
stage_bounds$celltype <- factor(stage_bounds$celltype, levels = stage_order)

pE <- ggplot() +
  geom_rect(
    data = stage_bounds,
    aes(xmin = pseudotime_q10, xmax = pseudotime_q90, ymin = -Inf, ymax = Inf, fill = celltype),
    alpha = 0.08, inherit.aes = FALSE
  ) +
  geom_ribbon(data = curve, aes(x = pseudotime_center, ymin = lo, ymax = hi), fill = "#80b1d3", alpha = 0.25) +
  geom_line(data = curve, aes(x = pseudotime_center, y = median), color = "#377eb8", linewidth = 0.7) +
  scale_fill_manual(values = stage_colors, drop = FALSE) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none") +
  labs(x = "Monocle pseudotime", y = "Local group-mixing entropy (0-1)")

# ---- Panel F: paired slopes ----
slopes <- read_csv(file.path(data_dir, "Fig2F_paired_slopes.csv"), show_col_types = FALSE)
slopes$group <- factor(slopes$group, levels = group_order)
slopes_long <- slopes %>%
  select(sample_id, group, entropy_early, entropy_late) %>%
  pivot_longer(cols = c(entropy_early, entropy_late), names_to = "stage", values_to = "entropy") %>%
  mutate(stage = recode(stage, entropy_early = "Early SPCs", entropy_late = "Late SPCs"))

pF <- ggplot(slopes_long, aes(x = stage, y = entropy, group = sample_id, color = group)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.2) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right") +
  labs(x = NULL, y = "Mean local group-mixing entropy (0-1)")

# ---- Panel G/H: GSEA PNGs (use existing rendered plots) ----
local_gsea_root <- file.path(data_dir, "gsea_pngs")
g_ks <- file.path(local_gsea_root, "Fig2G_CtrlvsKS_HALLMARK_OXIDATIVE_PHOSPHORYLATION.png")
g_az <- file.path(local_gsea_root, "Fig2G_CtrlvsAZFc_Del_HALLMARK_OXIDATIVE_PHOSPHORYLATION.png")
h_ks <- file.path(local_gsea_root, "Fig2H_CtrlvsKS_HALLMARK_APOPTOSIS.png")
h_az <- file.path(local_gsea_root, "Fig2H_CtrlvsAZFc_Del_HALLMARK_APOPTOSIS.png")

g_img <- function(path) {
  if (!file.exists(path)) {
    return(ggplot() + theme_void() + ggtitle(basename(path)))
  }
  if (requireNamespace("magick", quietly = TRUE)) {
    return(ggdraw() + draw_image(path))
  }
  if (requireNamespace("png", quietly = TRUE)) {
    img <- png::readPNG(path)
    grob <- grid::rasterGrob(img, interpolate = TRUE)
    return(ggdraw() + draw_grob(grob))
  }
  ggdraw() + draw_label("Install magick or png to render GSEA images")
}

pG <- g_img(g_ks) + g_img(g_az) + plot_layout(nrow = 1)
pH <- g_img(h_ks) + g_img(h_az) + plot_layout(nrow = 1)

# ---- save individual panels ----
ggsave(file.path(out_dir, "Fig2A.pdf"), pA, width = 8, height = 3.2)
ggsave(file.path(out_dir, "Fig2B.pdf"), pB, width = 4, height = 3.2)
ggsave(file.path(out_dir, "Fig2C.pdf"), pC, width = 8, height = 3.5)
ggsave(file.path(out_dir, "Fig2D.pdf"), pD, width = 10, height = 2.2)
ggsave(file.path(out_dir, "Fig2E.pdf"), pE, width = 6, height = 2.6)
ggsave(file.path(out_dir, "Fig2F.pdf"), pF, width = 4, height = 2.6)
ggsave(file.path(out_dir, "Fig2G.pdf"), pG, width = 6, height = 2.2)
ggsave(file.path(out_dir, "Fig2H.pdf"), pH, width = 6, height = 2.2)

# ---- assemble final figure ----
panelA <- tag_panel(pA, "A")
panelB <- tag_panel(pB, "B")
panelC <- tag_panel(pC, "C")
panelD <- tag_panel(pD, "D")
panelE <- tag_panel(pE, "E")
panelF <- tag_panel(pF, "F")
panelG <- tag_panel(pG, "G")
panelH <- tag_panel(pH, "H")

final <- (panelA | panelB) / panelC / panelD / (panelE | panelF) / (panelG | panelH)

ggsave(file.path(out_dir, "Figure2_rebuild.pdf"), final, width = 12, height = 14)
ggsave(file.path(out_dir, "Figure2_rebuild.png"), final, width = 12, height = 14, dpi = 300)

message("Saved panels + combined figure to: ", out_dir)
