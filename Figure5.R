#=====================================================================================
#
#  Code chunk 1 加载包、路径和数据
#
#=====================================================================================

rm(list=ls())


library(Cairo)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(pathview)
library(ggplot2)
library(forcats)
library(org.Hs.eg.db)

setwd("F:/10_R/Data_Output/Panel_4")
load("F:/10_R/Data_Output/Panel_3/DEG_all.RData")

colSAMPLE <- read.csv(
  "F:/10_R/Data_Input/colSAMPLE_Clean_LI.csv",
  header=T,row.names=1)

SAMPLE <- read.csv(
  "F:/10_R/Data_Input/RawCounts_Rearranged_Clean_LI.csv",
  header=T,row.names=1)


#=====================================================================================
#
#  Code chunk 3 富集分析（UCEC vs ConEC）
#
#=====================================================================================

genelist <- res[,1:2]

ENTREZID <- bitr(
  rownames(res), 
  fromType="SYMBOL", 
  toType = "ENTREZID", 
  OrgDb="org.Hs.eg.db"
)

genelist <- res[which(ENTREZID$SYMBOL %in% rownames(res)),]
genelist <- cbind(ENTREZID,genelist[,2])
rownames(genelist) <- genelist$ENTREZID

genelist <- genelist[,3]
names(genelist) <- ENTREZID$ENTREZID
genelist <- sort(genelist,decreasing = TRUE)

##
gsegocc <- gseGO(geneList     = genelist,#根据LogFC排序后的基因列表
                 OrgDb        = org.Hs.eg.db,
                 ont          = "CC",#GO分析的模块
                 minGSSize    = 10,#最小基因集的基因数
                 maxGSSize    = 500,#最大基因集的基因数
                 pvalueCutoff = 0.05,#p值的阈值
                 verbose      = FALSE)#是否输出提示信息

gsegocc1 <- simplify(
  gsegocc,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
gsegocc1 <- setReadable(gsegocc1,org.Hs.eg.db)
a <- gsegocc1@result[order(abs(gsegocc1@result$NES),decreasing=T)[1:50],]
a <- a[c(1,2,3,5,9,10,11,15,17,20,23,28,30,31,41,44,45,50),]
a <- a[-c(4,7,9,12,15,17),]

##
CairoPDF(file = "RP_GO_CC.pdf",8.5,4)
ggplot(a,
       aes(NES, fct_reorder(Description, NES), fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),trans = "log10", 
                       guide=guide_colorbar(reverse=TRUE, order = 1)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("Cellular Components")
dev.off()










gsegobp <- gseGO(geneList     = genelist,#根据LogFC排序后的基因列表
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",#GO分析的模块
                 minGSSize    = 10,#最小基因集的基因数
                 maxGSSize    = 500,#最大基因集的基因数
                 pvalueCutoff = 0.05,#p值的阈值
                 verbose      = FALSE)#是否输出提示信息

gsegobp1 <- simplify(
  gsegobp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

a <- gsegobp1@result[order(abs(gsegobp1@result$NES),decreasing=T)[1:180],]
a <- a[c(
  "GO:0032196","GO:0099560","GO:0007606",
  "GO:0006613","GO:0006413","GO:0000956",
  "GO:0006958","GO:0098760","GO:0060333",
  "GO:0048002","GO:0140467","GO:0070671","GO:0002227","GO:0032536"
),] 

##
CairoPDF(file = "RP_GO_BP.pdf",8.5,4)
ggplot(a,
       aes(NES, fct_reorder(Description, NES), fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),trans = "log10", 
                       guide=guide_colorbar(reverse=TRUE, order = 1)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("Biological Processes")
dev.off()




gsegomf <- gseGO(geneList     = genelist,#根据LogFC排序后的基因列表
                 OrgDb        = org.Hs.eg.db,
                 ont          = "MF",#GO分析的模块
                 minGSSize    = 10,#最小基因集的基因数
                 maxGSSize    = 500,#最大基因集的基因数
                 pvalueCutoff = 0.05,#p值的阈值
                 verbose      = FALSE)#是否输出提示信息

gsegomf1 <- simplify(
  gsegomf,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

a <- gsegomf1@result[order(abs(gsegomf1@result$NES),decreasing=T)[1:27],]
a <- a[c(
  "GO:0030021","GO:0033691","GO:0005549","GO:0003735",
  "GO:0042605","GO:0023026","GO:0042288","GO:0034987",
  "GO:0042834","GO:0015078","GO:0031492","GO:0019843"
),]

##
CairoPDF(file = "RP_GO_MF.pdf",10.5,4)
ggplot(a,
       aes(NES, fct_reorder(Description, NES), fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),trans = "log10", 
                       guide=guide_colorbar(reverse=TRUE, order = 1)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("Molecular Functions")
dev.off()



gsekegg <- gseKEGG(
  geneList     = genelist,
  organism     = 'hsa',
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  pAdjustMethod = "BH"
)
a <- gsekegg@result[order(abs(gsekegg@result$NES),decreasing=T)[1:85],]
a <- a[c(
  "hsa03010","hsa04612","hsa00190","hsa05100","hsa05322","hsa05203",
  "hsa05330","hsa04940","hsa04141","hsa04740","hsa05321","hsa04080"
  
),]
a <- setReadable(gsekegg, org.Hs.eg.db,keyType="ENTREZID")


CairoPDF(file = "RP_KEGG.pdf",7.5,4)
ggplot(a, 
       aes(NES, fct_reorder(Description, NES), fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),trans = "log10", 
                       guide=guide_colorbar(reverse=TRUE, order = 1)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("KEGG")
dev.off()

pathview(
  gene.data = genelist, 
  pathway.id = 'hsa04612',
  species="hsa", 
)

pathview(
  gene.data = genelist, 
  pathway.id = 'hsa04740',
  species="hsa", 
)

save.image("Panel_4_kai.RData")



CairoPDF(file = "GSEKEGG_Olfactory transduction.pdf",12,9)
gseaplot2(gsekegg, geneSetID = 'hsa04740',
          base_size = 20,
          pvalue_table = T,  #是否展示P值
          rel_heights = c(5,2,2),
          title = "Olfactory Transduction",
          ES_geom = 'line'
) 
dev.off()

CairoPDF(file = "GSEKEGG_Antigen processing and presentation.pdf",12,9)
gseaplot2(gsekegg, geneSetID = 'hsa04612',
          base_size = 20,
          pvalue_table = T,  #是否展示P值
          rel_heights = c(5,2,2),
          title = "Antigen Processing and Presentation",
          ES_geom = 'line'
) 
dev.off()




