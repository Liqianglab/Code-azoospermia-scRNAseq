
#ClusterGVis进行个性化热图可视化#
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
BiocManager::install("ComplexHeatmap")
BiocManager::install("clusterProfiler")
BiocManager::install("TCseq")
BiocManager::install("Mfuzz")
BiocManager::install("monocle3")
BiocManager::install("org.Mm.eg.db")
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages("circlize")
install.packages("Seurat")
install.packages("SeuratObject")
install.packages("Matrix")
library(SeuratObject)
library(Seurat)
library(monocle3)
library(circlize)
library(Matrix)


BiocManager::install("batchelor")
install.packages("batchelor")
install.packages("monocle")

options(timeout = 3000)  # 增加超时限制（单位为秒）

library(ClusterGVis)
# Note: please update your ComplexHeatmap to the latest version!
# install.packages("devtools")
devtools::install_github("junjunlab/ClusterGVis")

install.packages("Matrix", type = "binary")
install.packages("spam", type = "binary")
install.packages("spdep", type = "binary")
install.packages("igraph", type = "binary")

################################################################################################
################################################################################################

Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("/Users/xq/Desktop/睾丸文章总/最新内容汇总/figure汇总/figure2/reactome/未命名文件夹")

# 读取 CSV 文件
data <- read.csv("/Users/xq/Desktop/睾丸文章总/最新内容汇总/figure汇总/figure2/reactome/未命名文件夹/average_expression_matrix_ST.csv", header = T, stringsAsFactors = FALSE)

# 去掉重复的第一列行
data <- data[!duplicated(data[, 1]), ]
rownames(data) <- data[, 1]
data <- data[, -1] 
# 显示数据框的前几行，检查结果
head(data)

# 显示数据的前几行
head(data)

library(SeuratObject)
library(Seurat)
library(monocle3)
library(circlize)
library(Matrix)
library(ClusterGVis)

class(data)
library(Biobase)

# 将你的数据转换为 ExpressionSet
data1 <- ExpressionSet(assayData = as.matrix(data))  # 假设 'data' 是你的数据框或矩阵
library(Biobase)
library(Mfuzz)

class(data1)
head(data1)
detach(package:monocle3, unload = TRUE)
exprs_data <- Biobase::exprs(data1)
# 执行 mfuzz 聚类
cm <- clusterData(exprs_data,
                  cluster.method = "mfuzz",
                  cluster.num = 2)

# using mfuzz for clustering
cm <- clusterData(exprs_data,
                  cluster.method = "mfuzz",
                  cluster.num = 6)

# using TCseq for clustering
ct <- clusterData(exprs_data,
                  cluster.method = "TCseq",
                  cluster.num = 6)

# using kemans for clustering
ck <- clusterData(exprs_data,
                  cluster.method = "kmeans",
                  cluster.num = 2)
# plot line only
visCluster(object = cm,
           plot.type = "line")
#可视热图
# plot heatmap only
visCluster(object = ck,
           plot.type = "heatmap")

# 添加感兴趣的基因 add gene name
markGenes = rownames(exprs_data)[sample(1:nrow(exprs_data),30,replace = F)]

pdf('addgene.pdf',height = 10,width = 6,onefile = F)
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes)
dev.off()
#添加趋势线
pdf('testbxcol.pdf',height = 10,width = 6)
visCluster(object = ck,
           plot.type = "both",
           column_names_rot = 45,
           add.box = T,
           add.line = F,
           boxcol = ggsci::pal_npg()(8))
dev.off()
# multiple groups
pdf('termlftsmp.pdf',height = 10,width = 11,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
      
           genes.gp = c('italic',fontsize = 12,col = "black"),
           annoTerm.data = termanno2,
           line.side = "left",
           go.col = rep(ggsci::pal_d3()(8),each = 3),
           go.size = "pval",
           mulGroup = c(2,2,2),
           mline.col = c(ggsci::pal_lancet()(3)))
dev.off()


# load GO term data
data("termanno2")

# check
head(termanno2,3)
#   id                               term     pval
# 1 C1              developmental process 3.17e-69
# 2 C1   anatomical structure development 1.44e-68
# 3 C1 multicellular organism development 1.36e-66

# adjust term colors and text size
pdf('termlfts.pdf',height = 10,width = 12,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           column_names_rot = 45,

           genes.gp = c('italic',fontsize = 12,col = "black"),
           annoTerm.data = termanno2,
           line.side = "left",
           go.col = rep(ggsci::pal_d3()(8),each = 3),
           go.size = "pval")
dev.off()



#同时添加GO和KEGG 结果
# using mfuzz for clustering
# mfuzz
cm <- clusterData(exprs_data,
                  cluster.method = "mfuzz",
                  cluster.num = 3)

library(org.Hs.eg.db)
# GO enrich for clusters
enrich.go <- enrichCluster(object = cm,
                           OrgDb = org.Hs.eg.db,
                           type = "BP",
                           pvalueCutoff = 0.05,
                           topn = 5,
                           seed = 5201314)


# check
head(enrich.go,3)

#            group                       Description       pvalue     ratio
# GO:0044282    C1  small molecule catabolic process 3.376035e-32 10.683761
# GO:0046395    C1 carboxylic acid catabolic process 8.123713e-20  6.623932
# GO:0016054    C1    organic acid catabolic process 1.057149e-19  6.623932

# KEGG enrich for clusters
enrich.kegg <- enrichCluster(object = cm,
                             OrgDb = org.Hs.eg.db,
                             type = "KEGG",
                             organism = "hsa",
                             pvalueCutoff = 0.05,
                             topn = 15,
                             seed = 5201314)

# 假设 cm 是通过 mfuzz 聚类得到的结果



# check
head(enrich.kegg,3)

#              group                                      Description       pvalue    ratio
# mmu04146        C1          Peroxisome - Mus musculus (house mouse) 1.087835e-16 8.880309
# mmu01200...2    C1   Carbon metabolism - Mus musculus (house mouse) 3.447857e-14 9.266409
# mmu00620        C1 Pyruvate metabolism - Mus musculus (house mouse) 1.466749e-10 5.019305
# plot

pdf('gokegg.pdf',height = 10,width = 16,onefile = F)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
        
           annoTerm.data = enrich.go,
           go.col = rep(jjAnno::useMyCol("calm",n = 6),each = 15),
           annoKegg.data = enrich.kegg,
           kegg.col = rep(jjAnno::useMyCol("stallion",n = 6),each = 15),
           line.side = "left")
dev.off()
pdf('go.pdf',height = 10,width = 16,onefile = F)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           
           annoTerm.data = enrich.go,
           go.col = rep(jjAnno::useMyCol("calm",n = 2),each = 6),
           
           line.side = "left")
dev.off()
# add gene name
markGenes = rownames(exprs_data)[sample(1:nrow(exprs_data),30,replace = F)]

# plot
pdf('gokegg.pdf',height = 10,width = 16,onefile = F)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich.go,
           go.col = rep(jjAnno::useMyCol("calm",n = 6),each = 5),
           annoKegg.data = enrich.kegg,
           kegg.col = rep(jjAnno::useMyCol("stallion",n = 6),each = 5),
           line.side = "left",
           sample.group = rep(c("sample1","sample2","sample3"),each = 2))
dev.off()





cm <- clusterData(exprs_data,
                  cluster.method = "mfuzz",
                  cluster.num = 4)

library(org.Hs.eg.db)
# GO enrich for clusters
enrich.go <- enrichCluster(object = cm,
                           OrgDb = org.Hs.eg.db,
                           type = "BP",
                           pvalueCutoff = 0.05,
                           topn = 5,
                           seed = 5201314)

pdf('go.pdf',height = 10,width = 16,onefile = F)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           
           annoTerm.data = enrich.go,
           go.col = rep(jjAnno::useMyCol("calm",n = 4),each = 5),
           
           line.side = "left")
dev.off()
#添加柱状图#
pdf('bar.pdf',height = 10,width = 16,onefile = F)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = T,
           genes.gp = c('italic',fontsize = 12,col = "black"),
           annoTerm.data = enrich.go,
           go.col = rep(ggsci::pal_d3()(4),each = 5),
           go.size = "pval",
           add.bar = T)
dev.off()

