#=====================================================================================
#
#  Code chunk 1 加载包和路径
#
#=====================================================================================

rm(list=ls())
library(WGCNA)
library(Cairo)
enableWGCNAThreads(8)

setwd("F:/10_R/Data_Output/Clean_LI/WGCNA")


#=====================================================================================
#
#  Code chunk 2 加载表达矩阵和表型信息
#
#=====================================================================================

datExpr <- t(read.csv(
  "F:/10_R/Data_Input/NormData_Rearranged_Clean_LI.csv", 
  header=T,row.names=1))  #注意表达矩阵与DEG不同，需转置
datTraits <- read.csv(
  "F:/10_R/Data_Input/datTraits_Clean_LI_bygroup.csv", 
  header=T,row.names=1)

nSamples = nrow(datExpr)  #样本数目
nGenes = ncol(datExpr)  #基因数目


#=====================================================================================
#
#  Code chunk 3 数据过滤（可选）
#
#=====================================================================================

## 选取MAD最大的前5000基因，减少运算量，但会丢失信息

datExpr = datExpr[,
                  order(
                    apply(datExpr,1,mad), decreasing = T
                  )
                  [1:5000]
]

collectGarbage()


#=====================================================================================
#
#  Code chunk 4 绘制样本聚类树与表型热图（可选）
#
#=====================================================================================

## 对样本聚类
sampleTree = hclust(dist(datExpr), method = "average")
## 赋予表型颜色，灰色代表NA
traitColors = numbers2colors(datTraits, signed = F)

## 绘图，使用Cairo保存
CairoPDF("Sample_dendrogram_trait_heatmap.pdf" , 8 , 6)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits), 
                    cex.colorLabels = 0.8,
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 5, 3, 1),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#=====================================================================================
#
#  Code chunk 5 计算软阈值
#
#=====================================================================================

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(
  datExpr, 
  powerVector = powers, 
  verbose = 5,
  networkType = "signed"  #官方教程建议signed
)


#=====================================================================================
#
#  Code chunk 6 生成软阈值散点图
#
#=====================================================================================

## power已存入sft$powerEstimate，水张图，无卵用
## 此图需要手动保存

sizeGrWindow(9, 5)
par(mfrow = c(1,2))

# Scale-free topology fit index as a function of the soft-thresholding power
plot(
  sft$fitIndices[,1], 
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[,1], 
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,cex=0.9,
  col="red"
)

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type="n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels=powers,
  cex=0.9,
  col="red"
)
abline(h=100,col="red")


#=====================================================================================
#
#  Code chunk 7 根据软阈值计算基因模块
#
#=====================================================================================

## 耗时40min，结果千万务必要保存！！！

#若已计算过此步骤，从下面的代码载入net
#load("net.RData")

net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 19000,
  networkType = "signed",
  TOMType = "signed", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)

## 千万务必要保存！！！
save(net, file = "net.RData")


#=====================================================================================
#
#  Code chunk 8 基因模块可视化
#
#=====================================================================================

## 重要但无卵用的水图，然而文献都喜欢放

# 生成与模块相对应的颜色
moduleColors = labels2colors(net$colors)

# 对模块及其聚类树绘图
CairoPDF("Cluster_Dendrogram.pdf" , 8 , 6)
plotDendroAndColors(
  net$dendrograms[[1]], 
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()

#=====================================================================================
#
#  Code chunk 9 保存必要结果
#
#=====================================================================================

rm("traitColors", "powers")
save.image("WGCNA_Episode_1.RData")


#=====================================================================================
#
#  The End of WGCNA Episode 1
#
#=====================================================================================