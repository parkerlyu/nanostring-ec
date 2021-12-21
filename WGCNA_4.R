#=====================================================================================
#
#  Code chunk 1 加载包、路径和数据
#
#=====================================================================================

rm(list=ls())
library(pheatmap)
library(clusterProfiler)
library(ReactomePA)
#library(DOSE)
library(Cairo)

setwd("F:/10_R/Data_Output/Clean_LI/WGCNA")

load("WGCNA_Episode_2.RData")


#=====================================================================================
#
#  Code chunk 2 指定感兴趣模块
#
#=====================================================================================

module.selected <- "green3"

gene.selected <- colnames(datExpr)[(moduleColors==module.selected)]


#=====================================================================================
#
#  Code chunk 3 绘制模块基因表达量热图
#
#=====================================================================================

## 将表达矩阵再次转置并取出指定的基因表达信息
dat.selected <- t(datExpr)[gene.selected,]

## 导入colSAMPLE并指定分组颜色
colSAMPLE <- read.csv(
  "F:/10_R/Data_Input/colSAMPLE_Clean_LI.csv",
  header=T,row.names=1)

annotation_col <- colSAMPLE[,c("cell","subject")]
ann_colors = list(colSAMPLE$cell, colSAMPLE$subject)

## 用pheatmap绘制热图
CairoPDF(paste("Heatmap_module_selected" , module.selected , ".pdf", sep = "") ,
  8 , 8)

pheatmap(
  dat.selected,
  scale = "row" , 
  cluster_cols = F,
  cluster_rows = T,
  clustering_distance_rows = "correlation",
  
  show_rownames=T,
  show_colnames=T,
  
  border = F,
  
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  
  treeheight_row = 15, 
  treeheight_col = 10,
  
  legend = T,
  
  annotation_col = annotation_col,
  annotation_legend = TRUE,
  annotation_colors = ann_colors,
  
  fontsize = 7,
  cellwidth = 15,
  cellheight = 6
) 

dev.off()


#=====================================================================================
#
#  Code chunk 4 GO富集分析
#
#=====================================================================================

## 将ENTREZID加入表达矩阵
ENTREZID <- bitr(
  rownames(dat.selected), 
  fromType="SYMBOL", 
  toType = "ENTREZID", 
  OrgDb="org.Hs.eg.db"
  )

## 利用基因ID进行GO富集分析，需联网，运行较慢
enrichBP <- enrichGO(
  ENTREZID$ENTREZID, 
  OrgDb = org.Hs.eg.db, 
  ont='BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2,
  keyType = 'ENTREZID',
  readable = T
  )

## 使用simplify必须要指定ontology
a <- simplify(enrichBP, cutoff=0.7, by="p.adjust", select_fun=min)

enrichCC <- enrichGO(
  ENTREZID$ENTREZID, 
  OrgDb = org.Hs.eg.db, 
  ont='CC',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2,
  keyType = 'ENTREZID',
  readable = T
)

## 使用simplify必须要指定ontology
b <- simplify(enrichCC, cutoff=0.7, by="p.adjust", select_fun=min)


library(ggplot2)
library(forcats)
library(DOSE)

ego3 <- mutate(enrichGO@result, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ego3 <- ego3[c(1,2,4,5,6,9,10,14,16,40,33,32),]

CairoPDF(file = "WGCNA_GO.pdf",8.5,4)
ggplot(ego3, showCategory = 12,
      aes(richFactor, fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(4, 9)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  ylab(NULL) +
  ggtitle("GO Enrichment Analysis")
dev.off()





#=====================================================================================
#
#  Code chunk 5 KEGG富集分析（太少了不能用）
#
#=====================================================================================

## 与GO富集分析相似
enrichKEGG <- enrichKEGG(
  ENTREZID$ENTREZID, 
  organism = 'hsa', 
  keyType = 'kegg', 
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH', 
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE,
#  readable = T
)

dotplot(enrichKEGG, showCategory=30)

cnetplot(
  enrichKEGG, 
  showCategory = 12, 
  categorySize="pvalue", 
  colorEdge = T, 
  circular = F
  )

heatplot(enrichKEGG)


#=====================================================================================
#
#  Code chunk 6 Reactome富集分析
#
#=====================================================================================

enrichPathway <- enrichPathway(
  ENTREZID$ENTREZID, 
  pvalueCutoff = 0.05, 
  readable = T
  )

## 此处可直接查看结果（其实是我还没弄好可视化方法）
enrichPathway@result


#=====================================================================================
#
#  Code chunk 7 DO富集分析
# 
#=====================================================================================

enrichDO <- enrichDO(
  ENTREZID$ENTREZID,
  ont = "DO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
#  universe      = names(geneList),
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = T
  )

## 可视化方法也没弄好，凑合看吧
enrichDO@result


#=====================================================================================
#
#  Code chunk 8 保存结果
# 
#=====================================================================================

save.image(
  paste("WGCNA_Result_", module.selected, ".RData")
)


#=====================================================================================
#
#  End of WGCNA part IV: The Final Episode
# 
#=====================================================================================