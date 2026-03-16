
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

setwd("/Users/xq/Desktop/1")  # 修改为你的实际路径

# 加载必要的包
library(ggplot2)
library(readr)
library(ggpubr)

# 读取数据
df <- read_csv("expression (1).csv")
colnames(df)[colnames(df) == "gname"] <- "cell_type"

# 提取所有基因列名（排除 cell 与 cell_type）
gene_cols <- colnames(df)[!colnames(df) %in% c("cell", "cell_type")]

# 固定分组顺序
cell_order <- c("Ctrl", "AZFc_Del", "iNOA_B", "iNOA_S", "KS")

# 定义颜色
group_colors <- c(
  "Ctrl" = "#00a23f",
  "AZFc_Del" = "#ff1f1d",
  "iNOA_B" = "#a763ac",
  "iNOA_S" = "#0067aa",
  "KS" = "#ff7f00"
)

# 创建输出目录
dir.create("violin_plots", showWarnings = FALSE)

# 定义仅对 Ctrl 做的比较组
ctrl_comparisons <- list(
  c("Ctrl", "AZFc_Del"),
  c("Ctrl", "iNOA_B"),
  c("Ctrl", "iNOA_S"),
  c("Ctrl", "KS")
)

# 开始绘图循环
for (gene in gene_cols) {
  plot_data <- data.frame(expression = df[[gene]], cell_type = df$cell_type)
  plot_data <- subset(plot_data, expression > 0)
  plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_order)
  
  p <- ggplot(plot_data, aes(x = cell_type, y = expression, fill = cell_type)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      comparisons = ctrl_comparisons,
      step.increase = 0.15  # 逐层抬高星号避免重叠
    ) +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 16) +
    labs(
      title = NULL,
      x = "Cell Type",
      y = gene
    ) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  
  # 保存PDF图像
  ggsave(filename = paste0("violin_plots/", gene, "_violin.pdf"), plot = p, width = 6, height = 5)
}


#########代谢途径评分的小提琴图##########
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

setwd("/Users/xq/Desktop/1")  # 修改为你的目录

library(ggplot2)
library(readr)
library(ggpubr)

# 读取评分数据
df <- read_csv("B3_LCs_P20073103.diff_PRO.h5ad - Gene set enrichment15_expression.csv")

# 标准化列名
colnames(df)[colnames(df) == "corticotropin releasing hormone pathway"] <- "score"
colnames(df)[colnames(df) == "Major cell types"] <- "cell_type"

# 去除0值（如果有）
df <- subset(df, score > 0)

# 固定细胞类型顺序
df$cell_type <- factor(df$cell_type, levels = c("LC1", "LC2", "LC3"))

# 指定颜色
group_colors <- c(
  "LC1" = "#0066AB",
  "LC2" = "#00A340",
  "LC3" = "#FF1F1C"
)

# 比较所有两两组合
comparisons <- list(
  c("LC1", "LC2"),
  c("LC1", "LC3"),
  c("LC2", "LC3")
)

# 绘图
p <- ggplot(df, aes(x = cell_type, y = score, fill = cell_type)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = comparisons,
    step.increase = 0.15
  ) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 16) +
  labs(
    title = NULL,
    x = "Cell Type",
    y = "corticotropin releasing hormone score"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# 保存图像
ggsave("corticotropin releasing hormone_score_violin.pdf", plot = p, width = 6, height = 5)





##############################################

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

setwd("/Users/xq/Desktop/1")  # 修改为你的实际路径

library(ggplot2)
library(readr)
library(ggpubr)

# 读取文件
df <- read_csv("B3_LCs_P20073103.diff_PRO.h5ad - Gene set enrichment15_expression.csv")

# 重命名列
colnames(df)[colnames(df) == "corticotropin releasing hormone.csv"] <- "score"
colnames(df)[colnames(df) == "gname"] <- "cell_type"

# 去除0值（如果有）
df <- subset(df, score > 0)

# 固定分组顺序
cell_order <- c("Ctrl", "AZFc_Del", "iNOA_B", "iNOA_S", "KS")
df$cell_type <- factor(df$cell_type, levels = cell_order)

# 定义颜色
group_colors <- c(
  "Ctrl" = "#00a23f",
  "AZFc_Del" = "#ff1f1d",
  "iNOA_B" = "#a763ac",
  "iNOA_S" = "#0067aa",
  "KS" = "#ff7f00"
)

# 定义 Ctrl 对比组
ctrl_comparisons <- list(
  c("Ctrl", "AZFc_Del"),
  c("Ctrl", "iNOA_B"),
  c("Ctrl", "iNOA_S"),
  c("Ctrl", "KS")
)

# 创建输出目录
dir.create("violin_plots", showWarnings = FALSE)

# 绘图
p <- ggplot(df, aes(x = cell_type, y = score, fill = cell_type)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = ctrl_comparisons,
    step.increase = 0.15
  ) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 16) +
  labs(
    title = NULL,
    x = "Group",
    y = "corticotropin releasing hormone Score"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# 保存PDF图像
ggsave("corticotropin releasing hormon_violin.pdf", plot = p, width = 6, height = 5)





Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

setwd("/Users/xq/Desktop/1")  # ← 修改为你的路径

library(ggplot2)
library(readr)
library(ggpubr)

# 读取文件
df <- read_csv("B3_LCs_P20073103.diff_PRO.h5ad - Gene set enrichment17_expression.csv")

# 重命名列
colnames(df)[colnames(df) == "Luteinizing hormone.csv"] <- "score"
colnames(df)[colnames(df) == "Major cell types"] <- "cell_type"

# 去除0分值
df <- subset(df, score > 0)

# 固定分组顺序
df$cell_type <- factor(df$cell_type, levels = c("LC1", "LC2", "LC3"))

# 配色
group_colors <- c(
  "LC1" = "#0066AB",
  "LC2" = "#00A340",
  "LC3" = "#FF1F1C"
)

# 比较
comparisons <- list(
  c("LC1", "LC2"),
  c("LC1", "LC3"),
  c("LC2", "LC3")
)

# 创建输出目录
dir.create("violin_plots", showWarnings = FALSE)

# 绘图
p <- ggplot(df, aes(x = cell_type, y = score, fill = cell_type)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = comparisons,
    step.increase = 0.15
  ) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 16) +
  labs(
    title = NULL,
    x = "Cell Type",
    y = "Luteinizing hormone Score"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# 保存PDF图像
ggsave("Luteinizing hormone_violin.pdf", plot = p, width = 6, height = 5)


