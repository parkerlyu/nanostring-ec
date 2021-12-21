#=====================================================================================
#
#  Code chunk 1 加载包、路径和数据
#
#=====================================================================================

rm(list=ls())

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
#library(enrichplot)
#library(ggridges)
#library(ReactomePA)
library(Cairo)
library(DOSE)




#=====================================================================================
#
#  Code chunk 2 差异基因分析1（Con EC vs Con Epi）
#
#=====================================================================================


sample <- SAMPLE[,1:7]
colsample <- colSAMPLE[1:7,1:2]
dds <- DESeqDataSetFromMatrix(
    sample, 
    colsample,
    design= ~ cell
)
dds <- DESeq(dds)
res0 <- results(dds, contrast = c("cell", "EC", "Epithelium"))

normSAMPLE <- counts(dds, normalized=TRUE)
normSAMPLE.mad <- apply(normSAMPLE, 1, mad)
normSAMPLE <- normSAMPLE[order(normSAMPLE.mad, decreasing=T), ]
write.table(normSAMPLE, file = "normSAMPLE.csv", row.names = T, col.names = T, quote = FALSE,sep=',')

rld <- rlog(dds, blind=FALSE)
logSAMPLE <- assay(rld)
logSAMPLE <- logSAMPLE[order(normSAMPLE.mad, decreasing=T), ]
write.table(logSAMPLE, file = "logSAMPLE.csv", row.names = T, col.names = T, quote = FALSE,sep=',')

## padj为NA的替换为1
res0$padj[is.na(res0$padj)] <- 1
res0 <- res0[order(res0$log2FoldChange),]

write.table(res0, file="ConEC2ConEpi.csv", row.names = T, col.names = T, quote = FALSE, sep=',')

## 提取差异表达基因，padj<0.1，Log2FC>0.5！！！
deseq2.sig0 <- data.frame(subset(
    res0, 
    padj < 0.05
    & abs(log2FoldChange) > 1
))

write.table(deseq2.sig0, file = "ConEC2ConEpi.sig.csv", row.names = T, col.names = T, quote = FALSE,sep=',') 
save(dds, res0, deseq2.sig0, file = "DEG.RData")

#=====================================================================================
#
#  Code chunk 3 绘制差异基因火山图1（Con EC vs Con Epi）
#
#=====================================================================================

res0 <- data.frame(res0)

res0$threshold <- 
    factor(
        ifelse(res0$padj < 0.05 & abs(res0$log2FoldChange) >= 1, 
               ifelse(res0$log2FoldChange>= 1 ,'Up','Down'),'No Significance'),
        levels=c('Up','Down','No Significance')
    )

deseq2.sig0$threshold <- 
    factor(
        ifelse(deseq2.sig0$padj < 0.05 & abs(deseq2.sig0$log2FoldChange) >= 1, 
               ifelse(deseq2.sig0$log2FoldChange>= 1 ,'Up','Down'),'No Significance'),
        levels=c('Up','Down','No Significance')
    )

deseq2.sig0 <- data.frame(deseq2.sig0)
# res0$padj <- res0$padj + 1*10^(-5)
# deseq2.sig0$padj <- deseq2.sig0$padj + 1*10^(-5)

CairoPDF(file = "VP.ConEC2ConEpi.pdf",3.8,3.8)
ggplot(
    res0,
    aes(
        x = log2FoldChange,
        y = -log10(padj),
        color = threshold
    ) 
)+ 
    geom_point(size = 0.6, pch = 20)+
    scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
    ylab('-log10padj')+
    xlab('log2FoldC')+
    theme_bw()+hange
theme(
    legend.position="none",
    legend.title = element_blank()
)+
    geom_vline(xintercept=-1, linetype=2, colour="gray30") +
    geom_vline(xintercept=1, linetype=2, colour="gray30") +
    geom_hline(yintercept=-log10(0.05), linetype=2, colour="gray30") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    ) +
    scale_x_continuous(limit = c(-5.8,5.8)) +
    
    geom_text_repel(
        data = deseq2.sig0,
        aes(
            x = log2FoldChange,
            y = -log10(padj),
            label = rownames(deseq2.sig0)),
        size = 2.5,
        segment.color = "black", 
        show.legend = T,
        min.segment.length = 0.5,
        max.overlaps= 7,
        force = 10
    )
dev.off()

#=====================================================================================
#
#  Code chunk 4 绘制差异基因热图（Con EC vs Con Epi）
#
#=====================================================================================

annotation_col <- colSAMPLE[,c("cell","subject")]

ann_colors = list(
    colSAMPLE$cell)


SAMPLE.sig <- normSAMPLE[rownames(deseq2.sig0),]

pheatmap(
    SAMPLE.sig,
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
    
    fontsize = 6,
    cellwidth = 15,
    cellheight = 6,
    
    filename = "HM.sig.pdf"
) 


#=====================================================================================
#
#  Code chunk 5 GSEA富集分析
#
#=====================================================================================

genelist <- res0[,1:2]

ENTREZID <- bitr(
    rownames(res0), 
    fromType="SYMBOL", 
    toType = "ENTREZID", 
    OrgDb="org.Hs.eg.db"
)

genelist <- res0[which(ENTREZID$SYMBOL %in% rownames(res0)),]
genelist1 <- cbind(ENTREZID,genelist[,2])
rownames(genelist1) <- genelist1$ENTREZID

genelist <- genelist1[,3]
names(genelist) <- ENTREZID$ENTREZID
genelist <- sort(genelist,decreasing = TRUE)


gsegocc <- gseGO(geneList     = genelist,#根据LogFC排序后的基因列表
                 OrgDb        = org.Hs.eg.db,
                 ont          = "CC",#GO分析的模块
                 minGSSize    = 10,#最小基因集的基因数
                 maxGSSize    = 500,#最大基因集的基因数
                 pvalueCutoff = 0.05,#p值的阈值
                 verbose      = T)#是否输出提示信息

gsegocc1 <- simplify(
    gsegocc,
    cutoff = 0.6,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

a <- gsegocc1@result[order(abs(gsegocc1@result$NES),decreasing=T)[1:27],]
a <- a[-c(1:3,5:7,16,17),]
a <- a[1:12,]

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
    cutoff = 0.6,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

a <- gsegobp1@result[order(abs(gsegobp1@result$NES),decreasing=T)[1:65],]
a <- a[c(4,7,8,9,14,17,20,22,27,29,32,33,35,36,51,59,60,61),]
a <- a[-c(2,5,12,16,17,7),]

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
    cutoff = 0.6,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

a <- gsegomf1@result[order(abs(gsegomf1@result$NES),decreasing=T)[1:15],]
a <- a[-c(8,10,15),]

##
CairoPDF(file = "RP_GO_MF.pdf",8.5,4)
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
    pvalueCutoff = 0.1,
    verbose      = FALSE,
    pAdjustMethod = "BH"
)
a <- gsekegg@result[order(abs(gsekegg@result$NES),decreasing=T)[1:18],]
a <- a[-c(2),]
a <- a[1:12,]

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














