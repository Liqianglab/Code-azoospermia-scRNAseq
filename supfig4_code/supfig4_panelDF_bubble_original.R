#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

input_csv <- if (length(args) >= 1) {
  args[[1]]
} else {
  "../supfig4_source_data/raw_support/SupFig4F_ST3-diseasevsctrl.csv"
}

output_pdf <- if (length(args) >= 2) {
  args[[2]]
} else {
  "../supfig4_source_data/raw_support/SupFig4F_ST3_疾病vsctrl-NES-10_replot.pdf"
}

plot_title <- if (length(args) >= 3) {
  args[[3]]
} else {
  "NES-Top Enriched Pathways in ST3 (Disease vs Control)"
}

df <- read_csv(input_csv, show_col_types = FALSE)
required_cols <- c("Description", "pvalue", "NES", "group")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns: ", paste(missing_cols, collapse = ", ")))
}

df$Description <- str_replace_all(df$Description, "^KEGG_", "")
df$Description <- str_to_title(str_replace_all(df$Description, "_", " "))

pval_range <- range(df$pvalue, na.rm = TRUE)
nes_range <- range(df$NES, na.rm = TRUE)

top_df <- df %>%
  group_by(group) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

p <- ggplot(top_df, aes(x = Description, y = group)) +
  geom_point(aes(color = pvalue, size = NES)) +
  scale_color_gradient(low = "#336699", high = "#66CC66", limits = pval_range) +
  scale_size_continuous(range = c(3, 12), limits = nes_range) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = "Term",
    y = "Group",
    title = plot_title,
    size = "NES",
    color = "p-value"
  ) +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))

ggsave(output_pdf, plot = p, width = 8, height = 6)
cat("Saved:", output_pdf, "\n")
