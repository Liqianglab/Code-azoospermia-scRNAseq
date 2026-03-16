#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

base_dir <- "/Users/xq/Desktop/睾丸文章总/最新内容汇总/figure汇总/figure7"
anno_path <- file.path(base_dir, "Annotation_Row.csv")
expr_path <- file.path(base_dir, "late阶段所有基因表达.csv")
oxphos_score_path <- file.path(base_dir, "B3_shengzhi_P20073103.diff_PRO.h5ad - Extraction-late - ODX_expression.csv")
out_dir <- file.path(base_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading annotation: ", anno_path)
anno_raw <- read_csv(anno_path, show_col_types = FALSE) %>%
  rename(pathway = Group_Row, gene = Names_Row) %>%
  mutate(
    pathway = case_when(
      pathway == "Glycosis" ~ "Glycolysis",
      TRUE ~ pathway
    ),
    .row_order = row_number()
  )

message("Reading expression: ", expr_path)
expr <- read_csv(expr_path, show_col_types = FALSE) %>%
  rename(gene = Gene, ctrl = CTRL, disease = Disease)

eps <- 1e-6
expr <- expr %>%
  mutate(
    log2fc = log2((disease + eps) / (ctrl + eps)),
    mean_expr = (ctrl + disease) / 2
  )

df <- anno_raw %>%
  inner_join(expr, by = "gene") %>%
  mutate(
    pathway = factor(pathway, levels = c("Glycolysis", "TCA cycle", "Oxidative phosphorylation"))
  ) %>%
  arrange(pathway, .row_order) %>%
  mutate(
    gene = factor(gene, levels = rev(unique(gene)))
  )

write_csv(
  df %>% select(pathway, gene, ctrl, disease, log2fc),
  file.path(out_dir, "late_stage_metabolism_log2fc_table.csv")
)

message("Plot: gene-level log2FC heatmap")
heatmap_plot <- ggplot(df, aes(x = "Disease vs Ctrl", y = gene, fill = log2fc)) +
  geom_tile(width = 0.9, height = 0.9, color = "white", linewidth = 0.15) +
  facet_grid(pathway ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "#f7f7f7",
    high = "#d73027",
    midpoint = 0,
    name = expression(log[2] * "FC")
  ) +
  labs(x = NULL, y = NULL, title = "Late primary SPCs: metabolic genes (Disease vs Ctrl)") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 6),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    plot.title = element_text(face = "bold", hjust = 0)
  )

ggsave(
  filename = file.path(out_dir, "late_stage_metabolism_log2fc_heatmap.pdf"),
  plot = heatmap_plot,
  width = 7,
  height = 14,
  useDingbats = FALSE
)
ggsave(
  filename = file.path(out_dir, "late_stage_metabolism_log2fc_heatmap.png"),
  plot = heatmap_plot,
  width = 7,
  height = 14,
  dpi = 300
)

message("Plot: OXPHOS complex-level summary (log2FC mean across genes)")
oxphos <- df %>%
  filter(pathway == "Oxidative phosphorylation") %>%
  mutate(
    complex = case_when(
      str_detect(as.character(gene), "^MT-ND") ~ "Complex I",
      str_detect(as.character(gene), "^NDUF") ~ "Complex I",
      str_detect(as.character(gene), "^SDH") ~ "Complex II",
      str_detect(as.character(gene), "^UQCR") ~ "Complex III",
      str_detect(as.character(gene), "^MT-CYB$") ~ "Complex III",
      str_detect(as.character(gene), "^COX") ~ "Complex IV",
      str_detect(as.character(gene), "^MT-CO") ~ "Complex IV",
      str_detect(as.character(gene), "^ATP5") ~ "Complex V",
      str_detect(as.character(gene), "^MT-ATP") ~ "Complex V",
      str_detect(as.character(gene), "^COQ") ~ "CoQ biosynthesis",
      TRUE ~ "Other"
    )
  )

complex_summary <- oxphos %>%
  group_by(complex) %>%
  summarise(
    mean_log2fc = mean(log2fc, na.rm = TRUE),
    median_log2fc = median(log2fc, na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  mutate(
    complex = factor(
      complex,
      levels = c("Complex I", "Complex II", "Complex III", "Complex IV", "Complex V", "CoQ biosynthesis", "Other")
    )
  )

complex_plot <- ggplot(complex_summary, aes(x = complex, y = mean_log2fc, fill = mean_log2fc)) +
  geom_col(width = 0.75, color = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "#f7f7f7",
    high = "#d73027",
    midpoint = 0,
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = expression("Mean " * log[2] * "FC"),
    title = "Late primary SPCs: oxidative phosphorylation (complex-level summary)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0)
  )

ggsave(
  filename = file.path(out_dir, "late_stage_oxphos_complex_log2fc_summary.pdf"),
  plot = complex_plot,
  width = 7,
  height = 3.5,
  useDingbats = FALSE
)
ggsave(
  filename = file.path(out_dir, "late_stage_oxphos_complex_log2fc_summary.png"),
  plot = complex_plot,
  width = 7,
  height = 3.5,
  dpi = 300
)

if (file.exists(oxphos_score_path)) {
  message("Plot: late-stage oxidative phosphorylation score violin by group")
  score_df <- read_csv(oxphos_score_path, show_col_types = FALSE) %>%
    rename(
      cell = 1,
      score = 2
    ) %>%
    mutate(
      gname = factor(gname, levels = c("Ctrl", "OA", "AZFc_Del", "iNOA_B", "KS"))
    )

  score_plot <- ggplot(score_df, aes(x = gname, y = score, fill = gname)) +
    geom_violin(trim = FALSE, linewidth = 0.25, color = "grey25") +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, color = "grey15") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      x = NULL,
      y = "Module score",
      title = "Late primary SPCs: oxidative phosphorylation score"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0)
    )

  ggsave(
    filename = file.path(out_dir, "late_stage_oxphos_score_violin.pdf"),
    plot = score_plot,
    width = 6,
    height = 3.5,
    useDingbats = FALSE
  )
  ggsave(
    filename = file.path(out_dir, "late_stage_oxphos_score_violin.png"),
    plot = score_plot,
    width = 6,
    height = 3.5,
    dpi = 300
  )
}

message("Done. Output dir: ", out_dir)

