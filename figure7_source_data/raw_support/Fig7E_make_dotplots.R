#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

viridis5 <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")

group_order <- c("Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS")

read_long <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)

  expr_cols <- names(df)[str_detect(names(df), " raw expression value$")]
  if (length(expr_cols) == 0) {
    stop("No `* raw expression value` columns found in: ", path)
  }

  df %>%
    select(any_of(c("cell", "gname")), all_of(expr_cols)) %>%
    pivot_longer(cols = all_of(expr_cols), names_to = "gene", values_to = "raw_expr") %>%
    mutate(
      gene = str_replace(gene, " raw expression value$", ""),
      raw_expr = as.numeric(raw_expr),
      gname = factor(gname, levels = group_order)
    ) %>%
    filter(!is.na(gname))
}

summarize_dot <- function(df_long) {
  df_long %>%
    group_by(gname, gene) %>%
    summarise(
      pct_expressed = mean(raw_expr > 0, na.rm = TRUE) * 100,
      mean_log1p = mean(log1p(raw_expr), na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
}

make_dotplot <- function(df_sum, title, gene_order = NULL) {
  if (is.null(gene_order)) {
    gene_order <- sort(unique(as.character(df_sum$gene)))
  }
  df_sum <- df_sum %>%
    mutate(
      gene = factor(gene, levels = gene_order),
      gname = factor(as.character(gname), levels = group_order)
    )

  pct_breaks <- c(10, 25, 50, 75)

  ggplot(df_sum, aes(x = gname, y = gene)) +
    geom_point(aes(size = pct_expressed, colour = mean_log1p), alpha = 0.95) +
    scale_size(
      range = c(0.8, 12),
      limits = c(0, 100),
      breaks = pct_breaks,
      name = "% expressed"
    ) +
    scale_colour_gradientn(colors = viridis5, name = "Mean log1p\nexpression") +
    labs(x = "Group", y = "Gene", title = title) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 24, 12, 12)
    )
}

save_plot <- function(plot_obj, out_prefix, width = 10, height = 3.6) {
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}

main <- function() {
  base_dir <- "/Users/xq/Desktop/睾丸补图/figure7"
  plots_dir <- file.path(base_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  inputs <- tibble::tibble(
    file = c(
      "ST -ldha-slc16a3.csv",
      "late-spc-slc16a1-ldhb.csv",
      "round-spc-slc16a1-ldhb.csv"
    ),
    out_prefix = c(
      "ST_lactate_export_dotplot",
      "Late_primary_SPCs_lactate_uptake_dotplot",
      "Round_spermatids_lactate_uptake_dotplot"
    ),
    title = c(
      "Sertoli cells: lactate export markers",
      "Late primary spermatocytes: lactate uptake/entry markers",
      "Round spermatids: lactate uptake/entry markers"
    ),
    gene_order = I(list(
      c("SLC16A3", "LDHA"),
      c("SLC16A1", "LDHB"),
      c("SLC16A1", "LDHB")
    ))
  )

  plots <- vector("list", nrow(inputs))

  for (i in seq_len(nrow(inputs))) {
    in_path <- file.path(base_dir, inputs$file[i])
    df_long <- read_long(in_path)
    df_sum <- summarize_dot(df_long)
    p <- make_dotplot(df_sum, inputs$title[i], gene_order = inputs$gene_order[[i]])

    out_prefix <- file.path(plots_dir, inputs$out_prefix[i])
    save_plot(p, out_prefix)
    plots[[i]] <- p
  }

  combined <- (plots[[1]] / plots[[2]] / plots[[3]]) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")

  save_plot(
    combined,
    file.path(plots_dir, "Figure7_lactate_shuttle_dotplots_combined"),
    width = 10,
    height = 11
  )
}

main()

