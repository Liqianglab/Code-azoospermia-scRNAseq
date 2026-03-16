# 设置工作目录
setwd("/Users/xq/Desktop/体细胞figure/44.大类亚群组间差异基因及功能富集分析/Diff/immune")

# 安装必要包（只需首次执行）
# install.packages(c("clusterProfiler", "enrichplot", "ggplot2"))
# BiocManager::install("msigdbr")

# 加载包
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)

# 读取数据
df <- read.csv("immue—INOASvsCTRL.csv")

# 构建 Log2FC 排序向量，必须是 named vector，按降序排列
gene_list <- df$Log2FC
names(gene_list) <- df$Gene.id
gene_list <- sort(gene_list, decreasing = TRUE)

# 获取 hallmark gene sets（MSigDB H）
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# 执行 GSEA 分析
gsea_result <- GSEA(gene_list,
                    TERM2GENE = hallmark[, c("gs_name", "gene_symbol")],
                    pvalueCutoff = 0.05)

# 绘制 dotplot（经典 GSEA 点图）
dotplot(gsea_result, showCategory = 20, font.size = 10, title = "Hallmark Pathways")

# 或者绘制 ridgeplot（分布图）
ridgeplot(gsea_result, showCategory = 10) + labs(title = "Enrichment Ridge Plot")
