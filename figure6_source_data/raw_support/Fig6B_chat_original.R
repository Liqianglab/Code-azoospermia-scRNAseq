
setwd("/Users/xq/Desktop/未命名文件夹")
library(tidygraph)
library(ggraph)
library(igraph)
library(readr)
library(dplyr)
library(gridExtra)

# ==== 读取数据 ====
df_ctrl <- read_csv("all-ctrl.csv") %>% mutate(group = "CTRL")
df_disease <- read_csv("all-disease.csv") %>% mutate(group = "DISEASE")

# ==== 合并数据，分别处理 ====
plot_group <- function(df, title_label, line_type = "solid") {
  df <- df[df$interaction_count > 0, ]
  df$source_color <- df$source
  nodes <- data.frame(name = unique(c(df$source, df$target)))
  graph <- tbl_graph(nodes = nodes, edges = df, directed = TRUE)
  
  node_order <- c("LCs_a", "LCs_b", "LCs_c", "STs")
  group.colors <- c(
    "LCs_a" = "#1F77B4",
    "LCs_b" = "#2CA02C",
    "LCs_c" = "#E74C3C",
    "STs"   = "#9467BD"
  )
  
  layout <- create_layout(graph, layout = "circle")
  node_pos <- layout %>% select(name, x, y)
  
  # ==== 计算中点并偏移 ====
  edges_label <- df %>%
    left_join(node_pos, by = c("source" = "name")) %>%
    rename(x_source = x, y_source = y) %>%
    left_join(node_pos, by = c("target" = "name")) %>%
    rename(x_target = x, y_target = y) %>%
    rowwise() %>%
    mutate(
      x_mid = (x_source + x_target) / 2,
      y_mid = (y_source + y_target) / 2,
      dx = x_target - x_source,
      dy = y_target - y_source,
      length = sqrt(dx^2 + dy^2),
      nx = -dy / length,
      ny = dx / length,
      x = x_mid + 0.2 * nx,
      y = y_mid + 0.2 * ny
    )
  
  # ==== 绘图 ====
  p <- ggraph(layout) +
    geom_edge_arc(
      aes(width = interaction_count, color = source_color),
      strength = 0.3,
      end_cap = circle(3, 'mm'),
      arrow = arrow(length = unit(3, 'mm'), type = "closed"),
      linetype = line_type,
      show.legend = FALSE
    ) +
    geom_text(
      data = edges_label,
      aes(x = x, y = y, label = interaction_count),
      size = 3.5, fontface = "bold", color = "black"
    ) +
    geom_node_point(aes(x = x, y = y, color = name), size = 8, show.legend = FALSE) +
    geom_node_text(aes(x = x, y = y, label = name), size = 4.5, fontface = "bold") +
    scale_color_manual(values = group.colors) +
    scale_edge_color_manual(values = group.colors) +
    scale_edge_width(range = c(0.5, 5)) +
    coord_fixed() +
    theme_void() +
    ggtitle(title_label)
  return(p)
}

# ==== 生成左右图并保存 ====
pdf("all_CTRL_vs_DISEASE_side_by_side.pdf", width = 14, height = 7)
grid.arrange(
  plot_group(df_ctrl, "CTRL", line_type = "solid"),
  plot_group(df_disease, "DISEASE", line_type = "dashed"),
  ncol = 2
)
dev.off()
