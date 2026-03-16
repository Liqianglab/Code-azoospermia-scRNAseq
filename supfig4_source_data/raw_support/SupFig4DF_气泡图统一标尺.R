
# 1. 设置工作目录和文件列表
setwd("/Users/xq/Desktop/睾丸文章数据/ST/疾病的GSEA")  # 请替换成你的实际目录路径
# 加载所需库
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(readr)

# 读取数据
df <- read_csv("ST3-diseasevsctrl.csv")

# 清洗 Description 名称
df$Description <- str_replace_all(df$Description, "^KEGG_", "")
df$Description <- str_to_title(str_replace_all(df$Description, "_", " "))

# 提取范围
pval_range <- range(df$pvalue, na.rm = TRUE)
nes_range <- range(df$NES, na.rm = TRUE)

# 每个group中选p值最小前10项
top_df <- df %>%
  group_by(group) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  ungroup()

# 绘图，size 显示 NES
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
    x = "Term", y = "Group",
    title = "NES-Top Enriched Pathways in ST3 (Disease vs Control)",
    size = "NES", color = "p-value"
  ) +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))

# 保存图形
ggsave("ST3_疾病vsctrl-NES-10.pdf", plot = p, width = 8, height = 6)
