
# 加载必要包
library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
setwd("/Users/xq/Desktop/LC-GO汇总")
# 读取数据

# 读取数据
df <- read_csv("GO_enrichment (7)的副本.csv")

# 添加 -log10(pvalue)
df <- df %>%
  mutate(logp = -log10(pvalue))

# 设置组别颜色（PDF一致）
group_colors <- c(
  "AZFc" = "#1f77b4",
  "iNOA_B" = "#d62728",
  "iNOA_S" = "#2ca02c",
  "KS" = "#9467bd"
)

# 获取所有组最大logp值，确保x轴一致
max_logp <- max(df$logp, na.rm = TRUE)

# 函数：绘制单个组的图
plot_group <- function(group_name) {
  df_group <- df %>%
    filter(ONTOLOGY == group_name) %>%
    slice_min(order_by = pvalue, n = 10)
  
  ggplot(df_group, aes(x = logp, y = reorder(Description, logp))) +
    geom_bar(stat = "identity", fill = group_colors[group_name]) +
    xlim(0, max_logp) +  # 统一横坐标范围
    labs(
      title = group_name,
      x = "-log10(p-value)",
      y = NULL
    ) +
    theme_classic(base_size = 14)
}

# 分别生成4个图
p1 <- plot_group("AZFc")
p2 <- plot_group("iNOA_B")
p3 <- plot_group("iNOA_S")
p4 <- plot_group("KS")

# 如果需要一起排版展示，可使用 patchwork 或 cowplot
# install.packages("patchwork")
library(patchwork)
(p1 / p2 / p3 / p4) + plot_layout(ncol = 1)
