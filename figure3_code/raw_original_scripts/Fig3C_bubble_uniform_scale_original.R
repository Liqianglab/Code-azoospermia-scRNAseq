library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

# 1. 设置工作目录和文件列表
setwd("/Users/xq/Desktop/体细胞+GSEA+热图/气泡图")  # 请替换成你的实际目录路径
file_list <- list.files(pattern = "^table_PKD1_down_.*\\.csv$")
cell_names <- gsub("table_PKD1_down_(.*)\\.csv", "\\1", file_list)

# 2. 合并所有数据以统一图例范围
all_data <- lapply(seq_along(file_list), function(i) {
  df <- read.csv(file_list[i])
  df$celltype <- cell_names[i]
  return(df)
})
merged_df <- bind_rows(all_data)

# 3. 提取全局pvalue和setSize范围
pval_range <- range(merged_df$pvalue, na.rm = TRUE)
size_range <- range(merged_df$setSize, na.rm = TRUE)

# 4. 清洗term名称
clean_description <- function(desc) {
  desc <- str_replace_all(desc, "^HALLMARK_", "")
  str_to_title(str_replace_all(desc, "_", " "))
}

# 5. 每个细胞类型绘图
for (i in seq_along(file_list)) {
  df <- read.csv(file_list[i])
  df$Description <- clean_description(df$Description)
  df$celltype <- cell_names[i]
  
  # 每组中选pvalue最大前5项
  top_df <- df %>%
    group_by(group) %>%
    arrange(desc(pvalue)) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  # 图：横轴是term，纵轴是group
  p <- ggplot(top_df, aes(x = Description, y = group)) +
    geom_point(aes(color = pvalue, size = setSize)) +
    scale_color_gradient(low = "#336699", high = "#66CC66", limits = pval_range) +
    scale_size_continuous(range = c(3, 10), limits = size_range) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10)
    ) +
    labs(
      x = "Term", y = "Group",
      title = paste0("Top P-Value Terms per Group (Max): ", cell_names[i])
    ) +
    guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))
  
  # 保存PDF
  ggsave(filename = paste0("bubble_plot_", cell_names[i], "_pvalmax_rotated.pdf"),
         plot = p, width = 8, height = 4)
}
