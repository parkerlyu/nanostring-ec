#=====================================================================================
#
#  Code chunk 1 加载包、路径和数据
#
#=====================================================================================

rm(list=ls())
library(WGCNA)
library(Cairo)
enableWGCNAThreads(8)

setwd("F:/10_R/Data_Output/Clean_LI/WGCNA")

load("WGCNA_Episode_1.RData")


#=====================================================================================
#
#  Code chunk 2 全基因网络的可视化（慎用）
#
#=====================================================================================

## 运行时间极长，占用内存极大，不运行应直接跳过
load(net$TOMFiles, verbose=T)
TOM <- as.matrix(TOM)

dissTOM = 1-TOM
save(TOM, file = "TOM.RData")
rm(TOM)

plotTOM = dissTOM^7
save(dissTOM, file = "dissTOM.RData")
rm(dissTOM)

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
save(plotTOM, file = "plotTOM.RData")
# Call the plot function
save(moduleColors, file = "modulaColors.RData")
#selectColors = moduleColors[nGenes]

TOMplot(plotTOM, 
        net$dendrograms, 
        moduleColors, 
        nThreads = 0,
        main = "Network heatmap plot, all genes")


#=====================================================================================
#
#  Code chunk 3 网络的可视化（随机抽取部分基因）
#
#=====================================================================================


load(net$TOMFiles, verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM

## 随机选取部分基因
nSelect = 2000
set.seed(20)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

plotDiss = selectTOM^7
diag(plotDiss) = NA

CairoPDF("Network_heatmap_selected_genes.pdf" , 8 , 6)
TOMplot(
  plotDiss, 
  selectTree, 
  selectColors, 
  col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
  main = "Network heatmap for selected genes")
dev.off()

#=====================================================================================
#
#  Code chunk 4 筛选数据用于网络导出（可选）
#
#=====================================================================================

nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]


#=====================================================================================
#
#  Code chunk 5 准备导出数据
#
#=====================================================================================

module.selected = "green3"

geneName = colnames(datExpr)
geneSelected = geneName[(moduleColors==module.selected)]
geneSelected = (moduleColors==module.selected)


modTOM = TOM[geneSelected, geneSelected]
geneSelected = geneName[(moduleColors==module.selected)]
dimnames(modTOM) = list(geneSelected, geneSelected)
## 模块对应的基因关系矩阵 


#=====================================================================================
#
#  Code chunk 6 导出网络用于VISANT
#
#=====================================================================================

vis <- exportNetworkToVisANT(
  modTOM,
  file = paste("VisANTInput-", module, ".txt", sep=""),
  weighted = TRUE,
  threshold = 0
  )


#=====================================================================================
#
#  Code chunk 7 导出网络用于Cytoscape
#
#=====================================================================================

cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("Kai-CytoscapeInput-edges-", paste(module.selected, collapse="-"), ".txt", sep=""),
  nodeFile = paste("Kai-CytoscapeInput-nodes-", paste(module.selected, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0,
  nodeNames = geneSelected, 
  nodeAttr = geneTraitSignificance[geneSelected,]
  )


#=====================================================================================
#
#  Code chunk 8 保存必要数据
#
#=====================================================================================

rm("TOM", "dissTOM")
save.image("WGCNA_Episode_2.RData")


#=====================================================================================
#
#  The End of WGCNA Episode 2
#
#=====================================================================================