###############################
# ① 环境初始化
###############################
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

###############################
# ② 加载 R 包
###############################
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)

###############################
# ③ 工作目录 & 数据
###############################
setwd("/Users/xq/Desktop/临床数据")          # ← 修改
df <- read_csv("Violin_Input_generated.csv") # ← 修改

###############################
# ④ 颜色 & 比较组合
###############################
group_colors <- c("Control" = "#0066AB",
                  "OA"      = "#00A340",
                  "NOA"     = "#ff1f1d")

my_comparisons <- list(
  c("Control", "NOA"),
  c("Control", "OA"),
  c("NOA",     "OA")
)

dir.create("violin_plots", showWarnings = FALSE)

###############################
# ⑤ 循环绘图
###############################
for (metric_name in unique(df$Metric)) {
  
  sub_df <- df %>%
    filter(Metric == metric_name) %>%
    mutate(Group = factor(Group, levels = c("Control", "OA", "NOA")))
  
  # 计算 y 轴上限 = 数据最大值 ×1.4（给统计线留头顶空间）
  y_limit <- max(sub_df$Value, na.rm = TRUE) * 1.4
  
  p <- ggplot(sub_df, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
    geom_boxplot(width = 0.2, color = "black", fill = "white", outlier.shape = NA) +
    stat_compare_means(
      comparisons   = my_comparisons,
      method        = "wilcox.test",
      label         = "p.signif",   # 星号或 ns
      hide.ns       = FALSE,        # 显示 ns
      step.increase = 0.15          # 自动分层：每层抬 15%
    ) +
    scale_fill_manual(values = group_colors) +
    coord_cartesian(ylim = c(NA, y_limit)) +  # 拉高 y 轴上限
    theme_classic(base_size = 16) +
    labs(title = metric_name, x = "Group", y = "Value") +
    theme(
      legend.position = "none",
      axis.title      = element_text(size = 16, color = "black"),
      axis.text       = element_text(size = 14, color = "black"),
      axis.line       = element_line(color = "black"),
      axis.ticks      = element_line(color = "black")
    )
  
  ggsave(paste0("violin_plots/", metric_name, "_violin.pdf"),
         plot = p, width = 6, height = 5)
}

message("OK！显著性括号已自动分层，全部 PDF 在 violin_plots/ 中。")
