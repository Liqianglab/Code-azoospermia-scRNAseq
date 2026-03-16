setwd("/Users/xq/Desktop/体细胞+GSEA+热图/气泡图/EC")
all=read.csv('table_PKD1_down .csv', header = TRUE, sep = ",", quote = "\"")
library(ggplot2)
library(stringr)

# 清洗 Description 列，将 HALLMARK_ 去除，并转换为首字母大写格式
all$Description <- str_replace_all(all$Description, "^HALLMARK_", "")
all$Description <- str_to_title(str_replace_all(all$Description, "_", " "))
pdf("bubble_plot_LC.pdf", width = 5, height = 8)  # 你可以根据需要调整宽度和高度

ggplot(all, aes(group, Description)) +
  geom_point(aes(color = pvalue, size = setSize)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_color_gradient(low = "#336699", high = "#66CC66") +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 1))

dev.off()

