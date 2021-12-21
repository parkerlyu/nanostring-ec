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
#  Code chunk 2 计算模块与表型的相关性
#
#=====================================================================================

## 计算ME（Module Eigengenes）
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs)

## 计算模块基因（实为Eigen Genes）与表型的相关性
moduleTraitCor = cor(MEs, datTraits , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

## 准备相关性与p值的数据集用于绘制热图
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

## 绘图，手动保存
par(mar = c(12, 20, 5, 5))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships")
)

## 清除无用变量
rm("textMatrix")


#=====================================================================================
#
#  Code chunk 3 筛选code chunk 2中的显著结果，重新绘制热图（可选）
#
#=====================================================================================

moduleTraitPvalue.sig <- moduleTraitPvalue[
  which(
    rowSums(moduleTraitPvalue[,] < 0.05) > 0
    ), 
  ]

MEs.sig <- MEs[,rownames(moduleTraitPvalue.sig)]
moduleTraitCor.sig  <- moduleTraitCor[rownames(moduleTraitPvalue.sig),]

textMatrix.sig <- paste(signif(moduleTraitCor.sig, 2), "\n(",
                   signif(as.matrix(moduleTraitPvalue.sig), 1), ")", sep = "")
dim(textMatrix.sig) <- dim(moduleTraitCor.sig)

par(mar = c(12, 20, 5, 5))
labeledHeatmap(Matrix = moduleTraitCor.sig,
               xLabels = names(datTraits),
               yLabels = names(MEs.sig),
               ySymbols = names(MEs.sig),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix.sig,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships Selected")
)


#=====================================================================================
#
#  Code chunk 4 绘制模块相关性热图与聚类树
#
#=====================================================================================

## 将表型加入ME矩阵，不知所谓的操作，手动保存
MET = orderMEs(cbind(MEs, datTraits))

## 总图
sizeGrWindow(5,7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, 
                      "Eigengene adjacency heatmap", 
                      marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, 
                      xLabelsAngle= 90
)

## 聚类树
sizeGrWindow(6,6)
par(cex = 1.0)
plotEigengeneNetworks(MET, 
                      "Eigengene dendrogram", 
                      marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE
)

## 热图
sizeGrWindow(10,10)
par(cex = 0)
plotEigengeneNetworks(MET, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(8,8,8,8),
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90
)


#=====================================================================================
#
#  Code chunk 5 计算基因在每个模块的参与度与p值
#
#=====================================================================================

modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

## 将颜色名加“MM”前缀
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")


#=====================================================================================
#
#  Code chunk 6 并计算基因在感兴趣表型的参与度与p值
#
#=====================================================================================

name.trait.selected <- "CellEC"

trait.selected <- as.data.frame(datTraits[,name.trait.selected])
names(trait.selected) <- name.trait.selected

geneTraitSignificance <- as.data.frame(cor(datExpr, trait.selected, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

## 将颜色名加“GS”前缀
names(geneTraitSignificance) <- paste("GS.", names(trait.selected), sep="")
names(GSPvalue) <- paste("p.GS.", names(trait.selected), sep="")


#=====================================================================================
#
#  Code chunk 7 选取感兴趣模块绘制相关性散点图
#
#=====================================================================================

module.selected <- "green3"

column <- match(module.selected, modNames)
moduleGenes = moduleColors==module.selected

CairoPDF(paste("Module membership vs. gene significance for ", module.selected), 6, 6)
verboseScatterplot(
  abs(
    geneModuleMembership[
      moduleGenes, 
      match(module.selected, modNames)
      ]
    ),
  abs(
    geneTraitSignificance[
      moduleGenes, 1
      ]
    ),
  xlab = paste("Module Membership in", module.selected, sep = " "),
  ylab = paste("Gene significance for", module.selected, sep = " "),
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black",
  abline = T, abline.color = 1, abline.lty = 2
)
dev.off()

#=====================================================================================
#
#  Code chunk 8 保存必要结果
#
#=====================================================================================

save(
  datExpr, datTraits, MEs, moduleColors,
  file = "WGCNA_Episode_3.RData"
  )


#=====================================================================================
#
#  The End of WGCNA Episode 3
#
#=====================================================================================