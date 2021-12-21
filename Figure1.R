#=====================================================================================
#
#  Code chunk 1 加载包，路径和文件
#
#=====================================================================================

library(ggplot2)
library(amap)
library(gplots)
library(heatmap.plus)
library(Cairo)
library(ggpubr)


rm(list=ls())
setwd("F:/10_R/Data_Output/Panel_1")

## 先加载通过QC的Clean数据，用于绘制样本间聚类图
## 由于回肠末端数据为离群值，且每组只有1个样本，故剔除之
SAMPLE <- read.csv(
  "F:/10_R/Data_Input/NormData_Rearranged_Clean_LI.csv",
  header=T,row.names=1)
colSAMPLE <- read.csv(
  "F:/10_R/Data_Input/colSAMPLE_Clean_LI.csv",
  header=T,row.names=1)
normSAMPLE <- read.csv(
  "F:/10_R/Data_Input/NormData_Rearranged_Clean_LI.csv",
  header=T,row.names=1)


#=====================================================================================
#
#  Code chunk 2 绘制样本间聚类热图
#
#=====================================================================================

pearson_cor <- as.matrix(cor(SAMPLE, method="pearson"))
hc <- hcluster(t(SAMPLE), method="pearson")

CairoPDF("Heatmap_sample_cluster_subject.pdf", 6 , 6)

heatmap.2(
  pearson_cor, 
  Rowv=as.dendrogram(hc), 
  symm=T, 
  ColSideColors = c("#00DAE0", "#E199FF")[factor(colSAMPLE[,1])],#subject
  trace="none",
  col=colorRampPalette(c("navy", "white", "firebrick3"))(50), 
  margins=c(5,5), 
  #  dendrogram = "column",
  #  offsetRow = -0.3,
  #  cexRow = 0.8
)

dev.off()

CairoPDF("Heatmap_sample_cluster_cell.pdf", 6 , 6)

heatmap.2(
  pearson_cor, 
  Rowv=as.dendrogram(hc), 
  symm=T, 
  ColSideColors = c("#96CA00", "#FF9289")[factor(colSAMPLE[,2])],#cell
  trace="none",
  col=colorRampPalette(c("navy", "white", "firebrick3"))(50), 
  margins=c(5,5), 
  #  dendrogram = "column",
  #  offsetRow = -0.3,
  #  cexRow = 0.8
)

dev.off()
CairoPDF("Heatmap_sample_cluster_location.pdf", 6 , 6)
heatmap.2(
  pearson_cor, 
  Rowv=as.dendrogram(hc), 
  symm=T, 
  ColSideColors = c("#82B7FF", "#D3BA00")[factor(colSAMPLE[,3])],#location
  trace="none",
  col=colorRampPalette(c("navy", "white", "firebrick3"))(50), 
  margins=c(5,5), 
  #  dendrogram = "column",
  #  offsetRow = -0.3,
  #  cexRow = 0.8
)

dev.off()
#=====================================================================================
#
#  Code chunk 3 PCA
#
#=====================================================================================

mypca <- prcomp(
  normSAMPLE,
  scale = F,
  center = F
)

mypca.data <- data.frame(
  colSAMPLE,
  mypca$rotation[,1:2]
)

CairoPDF(file = "PCA_kai2.pdf",6,5)

#par(pin = c(5,5))
ggscatter(
  mypca.data,
  x = "PC1",
  y = "PC2",
  shape = "cell",
  color = "subject",
  #    ellipse = TRUE,
  #   ellipse.level = 0.80,
  size = 4,
#  alpha = 0.7,
  palette = c("navy", "firebrick3"),
  ggtheme = theme_dose(13)
  #,
  #  label = rownames(mypca.data),
  #  repel = TRUE
)

dev.off()


#=====================================================================================
#
#  Code chunk 4 全基因热图
#
#=====================================================================================

annotation_col <- colSAMPLE[,c("cell","subject")]

ann_colors = list(
  colSAMPLE$cell,
  colSAMPLE$subject
)



CairoPDF(paste("Heatmap_gene_uncluster_" , yellow.page[n,1] , ".pdf", sep = "") , 8 , 6)

pheatmap(
  SAMPLE,
  scale = "row" , 
  cluster_cols = F,
  cluster_rows = T,
  clustering_distance_rows = "correlation",
  
  show_rownames=F,  #Always false here
  show_colnames=T,  #Whether to show ROI names
  
  border = F,
  
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  
  treeheight_row = 15, 
  treeheight_col = 10,
  
  legend = T,
  
  annotation_col = annotation_col,
  annotation_legend = TRUE,
  annotation_colors = ann_colors,
  
  fontsize = 7,
  cellwidth = 15
  
) 

dev.off()








