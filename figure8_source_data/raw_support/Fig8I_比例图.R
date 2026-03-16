################################################################################################################################
#######比例图
#比例图
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
col_vector=c("#FF8AB6","#00A23F" ,"#FF1F1D" ,"#FF7F00" ,"#0067AA" ,"#A763AC" ,"#B45B5D")
pB2_df <- read.csv("/Users/xq/Desktop/睾丸文章总/最新内容汇总/figure汇总/figure5/体细胞表达/2.csv", header = TRUE, stringsAsFactors = FALSE)

#调整样本及cluster的顺序为指定顺序
cluster_order <- c("TCells","cDCs","LCs","ECs","Macrophages","Myoids","STs")
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = cluster_order)
sample_order <-c("Ctrl","OA")
pB2_df$Sample <- factor(pB2_df$Sample, levels = sample_order)

pB4 <- ggplot(data = pB2_df, aes(y =Number, x = Sample, fill =Cluster)) +
  geom_bar(stat = "identity", width=0.9,position="fill", color = "black", size = 0.2)+
  scale_fill_manual(values=col_vector) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB4
# 打印图表
print(pB4)
